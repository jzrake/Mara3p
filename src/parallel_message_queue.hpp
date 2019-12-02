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
     * @brief      Send a non-blocking message to a recipient process. The
     *             message buffer is stolen, and becomes a part of the Request
     *             object, which is in turn maintained in this message queue's
     *             container of send requests. Completed requests are removed
     *             when poll() is called.
     *
     * @param[in]  message         The message to send
     * @param[in]  recipient_rank  The message recipient rank
     */
    void push(buffer_t&& message, int recipient_rank)
    {
        send_requests.push_back(mpi::comm_world().isend(std::move(message), recipient_rank));
    }




    /**
     * @brief      Poll the send requests, deleting any that have completed.
     *             Then issue a non-blocking probe command to see if any MPI has
     *             any messages waiting for us. Return a vector all waiting
     *             messages.
     *
     * @return     A std::vector of messages waiting to be recieved
     */
    std::vector<buffer_t> poll()
    {
        send_requests.erase(std::remove_if(
            send_requests.begin(),
            send_requests.end(), std::mem_fn(&mpi::Request::is_ready)), send_requests.end());

        auto result = std::vector<buffer_t>();

        while (true)
        {
            auto comm = mpi::comm_world();
            auto status = comm.iprobe();

            if (status.is_null())
            {
                break;
            }
            result.push_back(comm.recv(status.source(), status.tag()));
        }
        return result;
    }

private:
    std::vector<mpi::Request> send_requests;
};
