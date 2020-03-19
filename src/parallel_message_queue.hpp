/**
 ==============================================================================
 Copyright 2019, Jonathan Zrake

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 ==============================================================================
*/




#pragma once
#include <algorithm>
#include <future>
#include <set>
#include <vector>
#include "parallel_mpi.hpp"




//=============================================================================
namespace mara {
    class MessageQueue;
}




//=============================================================================
class mara::MessageQueue
{
public:

    //=========================================================================
    using buffer_t = mpi::buffer_t;




    /**
     * @brief      Return true if there are no unsatisfied outgoing send
     *             requests, and all async pushes are completed.
     *
     * @return     True or false
     */
    bool empty() const
    {
        return send_requests.empty() && async_pushes.empty();
    }




    /**
     * @brief      Send a non-blocking message to a recipient process. The
     *             message buffer is stolen, and becomes a part of the Request
     *             object, which is in turn maintained in this message queue's
     *             container of send requests. Completed requests are removed
     *             when poll() is called.
     *
     * @param[in]  message         The message to send
     * @param[in]  recipient_rank  The message recipient rank
     */
    void push(buffer_t message, int recipient_rank, int tag=mpi::any_tag)
    {
        send_requests.push_back(mpi::comm_world().isend(std::move(message), recipient_rank, tag));
    }




    /**
     * @brief      Convenience method to push a message to a number of unique
     *             recipients.
     *
     * @param[in]  message          The message to push
     * @param[in]  recipient_ranks  The ranks of the recipient processes
     * @param[in]  tag              The message tag (optional)
     */
    void push(buffer_t message, std::set<int> recipient_ranks, int tag=mpi::any_tag)
    {
        for (auto recipient : recipient_ranks)
        {
            push(message, recipient, tag);
        }
    }




    /**
     * @brief      Push an asynchronous message to a set of recipients. Async
     *             pushes are polled, and the completed ones sent, each time
     *             check_outgoing is called.
     *
     * @param[in]  message          The message future
     * @param[in]  recipient_ranks  The ranks of the recipient processes
     * @param[in]  tag              The message tag (optional)
     */
    void push(std::future<buffer_t> message, std::set<int> recipient_ranks, int tag=mpi::any_tag)
    {
        async_pushes.emplace_back(std::move(message), std::move(recipient_ranks), tag);
    }




    /**
     * @brief      Poll the send requests, deleting any that have completed.
     *             Then issue a non-blocking probe command to see if MPI has any
     *             messages waiting for us. Return a vector of all received
     *             messages.
     *
     * @return     A std::vector of messages that were received
     */
    std::vector<buffer_t> poll()
    {
        auto result = std::vector<buffer_t>();

        for (auto&& [tag, buffer] : poll_tags())
        {
            result.push_back(std::move(buffer));
        }
        return result;
    }




    /**
     * @brief      Same as poll, but returns a vector of pair of tag-buffer
     *             pairs.
     *
     * @return     A std::vector of tags and messages that were received
     */
    std::vector<std::pair<int, buffer_t>> poll_tags()
    {
        auto result = std::vector<std::pair<int, buffer_t>>();

        while (true)
        {
            auto comm = mpi::comm_world();
            auto status = comm.iprobe();

            if (status.is_null())
            {
                break;
            }
            result.emplace_back(status.tag(), comm.recv(status.source(), status.tag()));
        }
        return result;
    }




    /**
     * @brief      Poll the asynchronous pushes and MPI non-blocking send
     *             requests. For each async push that is ready, the message
     *             is extracted from the future (rendering it invalid) and
     *             sent to the recipient ranks.
     */
    void check_outgoing()
    {
        for (auto& [future_message, recipients, tag] : async_pushes)
        {
            if (future_message.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
            {
                push(future_message.get(), recipients, tag);
            }
        }

        auto ar = std::remove_if(async_pushes.begin(), async_pushes.end(), [] (const auto& p) { return ! std::get<0>(p).valid(); });
        auto sr = std::remove_if(send_requests.begin(), send_requests.end(), std::mem_fn(&mpi::Request::is_ready));

        async_pushes.erase(ar, async_pushes.end());
        send_requests.erase(sr, send_requests.end());
    }




    /**
     * @brief      Return a const reference to the vector of open send requests.
     *
     * @return     The open send requests
     */
    const std::vector<mpi::Request>& open_send_requests() const
    {
        return send_requests;
    }




private:
    using async_push_t = std::tuple<std::future<buffer_t>, std::set<int>, int>;

    std::vector<mpi::Request> send_requests;
    std::vector<async_push_t> async_pushes;
};
