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
#include <cstring>
#include <functional>
#include <numeric>
#include <vector>
#include <mpi.h>




// ============================================================================
namespace mpi
{
    class Session;
    class Communicator;
    class Request;
    class Status;

    using buffer_t = std::vector<char>;

    inline Communicator comm_world();
    constexpr int any_tag = MPI_ANY_TAG;
    constexpr int any_source = MPI_ANY_SOURCE;

    enum class operation
    {
        max,
        min,
        sum,
        prod,
        land,
        lor,
        band,
        bor,
        maxloc,
        minloc,
    };

    inline auto get_op(operation op)
    {
        switch (op)
        {
            case operation::max:    return MPI_MAX;
            case operation::min:    return MPI_MIN;
            case operation::sum:    return MPI_SUM;
            case operation::prod:   return MPI_PROD;
            case operation::land:   return MPI_LAND;
            case operation::lor:    return MPI_LOR;
            case operation::band:   return MPI_BAND;
            case operation::bor:    return MPI_BOR;
            case operation::maxloc: return MPI_MAXLOC;
            case operation::minloc: return MPI_MINLOC;
        }
        throw std::invalid_argument("mpi::get_op (unknown MPI operation)");
    }

    inline void throw_unless_valid_tag(int tag)
    {
        if (tag >= MPI_TAG_UB)
        {
            throw std::overflow_error("tag >= MPI_TAG_UB");
        }
    }
}




// ============================================================================
class mpi::Session
{
public:
    Session()
    {
        MPI_Init(0, nullptr);
    }
    ~Session()
    {
        MPI_Finalize();
    }
    Session(const Session& other) = delete;
};




// ============================================================================
/**
 * A thin RAII wrapper around the MPI_Request struct. This is a movable, but
 * non-copyable object. Keep it around on the stack in order to check on the
 * status of a non-blocking communication and retrieve the message content.
 * Requests are cancelled and deallocated (if necessary) when they go out of
 * scope.
 */
class mpi::Request
{
public:


    /**
     * Default constructor, creates a null request.
     */
    Request() {}


    /**
     * Request is a unique object, no copy's are permitted.
     */
    Request(const Request& other) = delete;


    /**
     * Move constructor. Steals ownership of the other.
     */
    Request(Request&& other)
    {
        buffer = std::move(other.buffer);
        request = other.request;
        other.request = MPI_REQUEST_NULL;
    }


    /**
     * Destructor. Cancels the request if one is pending. For this reason, the
     * request returned by non-blocking communications must be retained on the
     * stack somewhere in order for the operation not to be cancelled.
     */
    ~Request()
    {
        cancel();
    }


    /**
     * Copy assignment is not permitted.
     */
    Request& operator=(const Request& other) = delete;


    /**
     * Move assignment. Cancels the current request (if one is pending) and
     * steals ownership of the other.
     */
    Request& operator=(Request&& other)
    {
        cancel();
        buffer = std::move(other.buffer);
        request = other.request;
        other.request = MPI_REQUEST_NULL;
        return *this;
    }


    /**
     * Cancel this request and reset its state to null.
     */
    void cancel()
    {
        if (! is_null())
        {
            if (! is_ready())
            {
                MPI_Cancel(&request);
            }
            MPI_Request_free(&request);
        }
    }


    /**
     * Return true if this request is null.
     */
    bool is_null() const
    {
        return request == MPI_REQUEST_NULL;
    }


    /**
     * Check to see whether the request has completed. Like above, except
     * this const method will not free the request if it has completed.
     */
    bool is_ready() const
    {
        int flag;
        MPI_Request_get_status(request, &flag, MPI_STATUS_IGNORE);
        return flag;
    }


    /**
     * Block until the request is fulfilled. After this method returns, the
     * get() method can be called to retrieve the message content.
     */
    void wait()
    {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }


    /**
     * Block until the message is completed. Then deallocate the request, and
     * return the message content.
     */
    const buffer_t& get()
    {
        wait();
        return buffer;
    }


    /**
     * Return the message content formatted as the given data type.
     */
    template <typename T>
    T get()
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");

        if (buffer.size() != sizeof(T))
        {
            throw std::logic_error("mpi::Request (received message has wrong size for data type)");
        }
        wait();

        auto value = T();
        std::memcpy(&value, &buffer[0], sizeof(T));
        return value;
    }


private:
    // ========================================================================
    friend class Communicator;
    MPI_Request request = MPI_REQUEST_NULL;
    buffer_t buffer;
};




// ============================================================================
/**
    This is a simple wrapper class around the MPI_Status struct.
*/
class mpi::Status
{
public:


    /**
     * Default-constructed status has is_null() == true.
     */
    Status() {}


    /**
     * Check for null-status. Null status is returned e.g. when iprobe does
     * not find any incoming messages:
     *
     *              if (comm.iprobe().is_null()) { }
     *
     */
    bool is_null() const
    {
        return null;
    }


    /**
     * Return the number of bytes in the message described by this status.
     */
    int count()
    {
        if (is_null())
        {
            return 0;
        }

        int count;
        MPI_Get_count(&status, MPI_CHAR, &count);
        return count;
    }


    /**
     * Get the rank of the message's source.
     */
    int source() const
    {
        if (is_null())
        {
            return -1;
        }
        return status.MPI_SOURCE;
    }


    /**
     * Get the tag of the message.
     */
    int tag() const
    {
        if (is_null())
        {
            return -1;
        }
        return status.MPI_TAG;
    }


private:
    // ========================================================================
    friend class Communicator;
    Status(MPI_Status status) : status(status), null(false) {}
    MPI_Status status;
    bool null = true;
};




// ============================================================================
class mpi::Communicator
{
public:


    /**
     * Default constructor, gives you MPI_COMM_NULL.
     */
    Communicator()
    {
    }


    /**
     * Communicator is a unique object, no copy's are permitted..
     */
    Communicator(const Communicator& other) = delete;


    /**
     * Move constructor, sets the other comm back to null.
     */
    Communicator(Communicator&& other)
    {
        comm = other.comm;
        owned = other.owned;
        other.comm = MPI_COMM_NULL;
    }


    /**
     * Destructor, closes the communicator unless it was null.
     */
    ~Communicator()
    {
        close();
    }



    /**
     * Move assignment. Steals the other communicator.
     */
    Communicator& operator=(Communicator&& other)
    {
        close();

        comm = other.comm;
        other.comm = MPI_COMM_NULL;
        other.owned = true;
        return *this;
    }


    /**
     * Close the communicator if it wasn't null.
     */
    void close()
    {
        if (! is_null())
        {
            if (owned)
            {
                MPI_Comm_free(&comm);
            }
            comm = MPI_COMM_NULL;
        }
    }


    /**
     * Return true if the communicator is null.
     */
    bool is_null() const
    {
        return comm == MPI_COMM_NULL;
    }


    /**
     * Return the number of ranks in the communicator. This returns zero for a
     * null communicator (whereas I think MPI implementations typically
     * produce an error).
     */
    int size() const
    {
        if (is_null())
        {
            return 0;
        }

        int res;
        MPI_Comm_size(comm, &res);
        return res;
    }


    /**
     * Return the rank of the communicator. This returns -1 for a null
     * communicator (whereas I think MPI implementations typically produce an
     * error).
     */
    int rank() const
    {
        if (is_null())
        {
            return -1;
        }

        int res;
        MPI_Comm_rank(comm, &res);
        return res;
    }


    /**
     * Block all ranks in the communicator at this points.
     */
    void barrier() const
    {
        MPI_Barrier(comm);
    }


    /**
     * Probe for an incoming message and return its status. This method blocks
     * until there is an incoming message to probe.
     */
    Status probe(int rank=any_source, int tag=any_tag) const
    {
        throw_unless_valid_tag(tag);
        MPI_Status status;
        MPI_Probe(rank, tag, comm, &status);
        return status;
    }


    /**
     * Probe for an incoming message and return its status. This method will
     * not block, but returns a null status if there was no message to probe.
     */
    Status iprobe(int rank=any_source, int tag=any_tag) const
    {
        throw_unless_valid_tag(tag);
        MPI_Status status;
        int flag;
        MPI_Iprobe(rank, tag, comm, &flag, &status);

        if (! flag)
        {
            return Status();
        }
        return status;
    }


    /**
     * Blocking-receive a message with the given source and tag. Return the
     * data as a buffer.
     */
    buffer_t recv(int source=any_source, int tag=any_tag) const
    {
        throw_unless_valid_tag(tag);
        auto status = probe(source, tag);
        auto buf = buffer_t(status.count(), 0);

        MPI_Recv(&buf[0], buf.size(), MPI_CHAR, source, tag, comm, MPI_STATUS_IGNORE);
        return buf;
    }


    /**
     * Non-blocking receive a message with the given source and tag. Return a
     * request object that can be queried for the completion of the receive
     * operation. Note that the request is cancelled if allowed to go out of
     * scope. You should keep the request somewhere, and call test(), wait(),
     * or get() in a little while.
     */
    Request irecv(int source=any_source, int tag=any_tag) const
    {
        throw_unless_valid_tag(tag);
        auto status = iprobe(source, tag);

        if (status.is_null())
        {
            return Request();
        }
        auto buf = buffer_t(status.count(), 0);

        MPI_Request request;
        MPI_Irecv(&buf[0], buf.size(), MPI_CHAR, source, tag, comm, &request);

        Request res;
        res.buffer = std::move(buf);
        res.request = request;
        return res;
    }


    /**
     * Blocking-send a buffer to the given rank.
     */
    void send(buffer_t buf, int rank, int tag=0) const
    {
        throw_unless_valid_tag(tag);
        MPI_Send(&buf[0], buf.size(), MPI_CHAR, rank, tag, comm);
    }


    /**
     * Blocking-sendrecv a buffer from dest to source
     */
    buffer_t sendrecv(buffer_t sendbuf, int dest, int source, int tag=0)
    {
        throw_unless_valid_tag(tag);
        auto status  = probe(source, tag);
        auto recvbuf = buffer_t(status.count(), 0);

        MPI_Sendrecv(&sendbuf[0], sendbuf.size(), MPI_CHAR, dest, tag,
                     &recvbuf, recvbuf.size(), MPI_CHAR, source, tag,
                     comm, MPI_STATUS_IGNORE);

        return recvbuf;
    }


    /**
     * Non-blocking send a buffer to the given rank. Returns a request object
     * that can be tested for completion or waited on. Note that the request
     * is cancelled if allowed to go out of scope. Also keep in mind your MPI
     * implementation may have chosen to buffer your message internally, in
     * which case the request will have completed immediately, and the
     * cancellation will have no effect. Therefore it is advisable to keep the
     * returned Request object, or at least do something equivalent to:
     *
     *              auto result = comm.isend("Message!", 0).get();
     *
     * Of course this would literally be a blocking send, but you get the
     * idea. In practice you'll probably store the request somewhere and check
     * on it after a while.
     */
    Request isend(buffer_t buf, int rank, int tag=0) const
    {
        throw_unless_valid_tag(tag);
        MPI_Request request;
        MPI_Isend(&buf[0], buf.size(), MPI_CHAR, rank, tag, comm, &request);

        Request res;
        res.buffer = std::move(buf);
        res.request = request;
        return res;
    }


    /**
     * Template version of a blocking send. You can pass any standard-layout
     * data type here.
     */
    template <typename T>
    void send(const T& value, int rank, int tag=0) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");
        auto buf = buffer_t(sizeof(T), 0);
        std::memcpy(&buf[0], &value, sizeof(T));
        send(buf, rank, tag);
    }


    /**
     * Template version of a non-blocking send. You can pass any standard-layout
     * data type here.
     */
    template <typename T>
    Request isend(const T& value, int rank, int tag=0) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");
        auto buf = buffer_t(sizeof(T), 0);
        std::memcpy(&buf[0], &value, sizeof(T));
        return isend(buf, rank, tag);
    }


    /**
     * Template version of a blocking receive. You can pass any
     * standard-layout data type here.
     */
    template <typename T>
    T recv(int rank, int tag=0) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");
        auto buf = recv(rank, tag);

        if (buf.size() != sizeof(T))
        {
            throw std::logic_error("received message has wrong size for data type");
        }

        auto value = T();
        std::memcpy(&value, &buf[0], sizeof(T));
        return value;
    }


    /**
     * Template version of a send-receive using blocking sends and
     * receives.
     *
     */
    template <typename T>
    T sendrecv(const T &send_value, int send_rank, int recv_rank, int tag=0) const
    {
        static_assert(std::is_trivially_copyable<T>::send_value, "type is not trivially copyable");

        auto sendbuf = buffer_t(sizeof(T));
        std::memcpy(&sendbuf[0], &send_value, sizeof(T));

        auto recvbuf = sendrecv(sendbuf, send_rank, recv_rank, tag);

        if (recvbuf.size() != sizeof(T))
        {
            throw std::logic_error("received message has wrong size for data type");
        }

        auto recv_value = T();
        std::memcpy(&recv_value, &recvbuf[0], sizeof(T));
        return recv_value;
    }


    /**
     * Execute a bcast operation with the given rank as the root.
     */
    template <typename T>
    void bcast(int root, T& value) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");
        MPI_Bcast(&value, sizeof(T), MPI_CHAR, root, comm);
    }


    /**
     * Execute a scatter communication with the given rank as root. The i-th
     * index of the send buffer is received by the i-th rank. The send buffer
     * is ignored by all processes except the root.
     */
    template <typename T>
    T scatter(int root, const std::vector<T>& values) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");

        if (root == rank() && values.size() != size())
        {
            throw std::invalid_argument("scatter send buffer must equal the comm size");
        }

        auto value = T();

        MPI_Scatter(
            &values[0], sizeof(T), MPI_CHAR,
            &value, sizeof(T), MPI_CHAR, root, comm);

        return value;
    }


    /**
     * Execute a scatter-v communication with the given rank as root. The i-th
     * index of the send buffer is received by the i-th rank. The send buffer
     * is ignored by all processes except the root.
     */
    template <typename T>
    std::vector<T> scatter(int root, const std::vector<std::vector<T>>& values) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");

        if (root == rank() && values.size() != size())
        {
            throw std::invalid_argument("scatter send buffer must equal the comm size");
        }

        if (rank() == root)
        {
            auto sendcounts = std::vector<int>(size());
            auto senddispls = std::vector<int>{0};
            auto sendbuf    = std::vector<T>();

            for (int i = 0; i < values.size(); ++i)
            {
                sendcounts[i] = values[i].size() * sizeof(T);
                sendbuf.insert(sendbuf.end(), values[i].begin(), values[i].end());
            }
            std::partial_sum(sendcounts.begin(), sendcounts.end() - 1, std::back_inserter(senddispls));

            auto recvcount = scatter(root, sendcounts);
            auto recvbuf   = std::vector<T>(recvcount / sizeof(T));

            MPI_Scatterv(
                &sendbuf[0], &sendcounts[0], &senddispls[0], MPI_CHAR,
                &recvbuf[0], recvcount, MPI_CHAR, root, comm);

            return recvbuf;
        }
        else
        {
            auto recvcount = scatter(root, std::vector<int>());
            auto recvbuf   = std::vector<T>(recvcount / sizeof(T));

            MPI_Scatterv(
                nullptr, nullptr, nullptr, MPI_CHAR,
                &recvbuf[0], recvcount, MPI_CHAR, root, comm);

            return recvbuf;
        }
    }


    /**
     * Execute an all-to-all communication with container of data. Each rank
     * sends the value at index i to rank i. The return value at index j
     * contains the character received from rank j.
     */
    template <typename T>
    std::vector<T> all_to_all(const std::vector<T>& sendbuf) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");

        if (sendbuf.size() != size())
        {
            throw std::invalid_argument("all_to_all send buffer must equal the comm size");
        }

        auto recvbuf = std::vector<T>(sendbuf.size(), T());

        MPI_Alltoall(
            &sendbuf[0], sendbuf.size() / size() * sizeof(T), MPI_CHAR,
            &recvbuf[0], recvbuf.size() / size() * sizeof(T), MPI_CHAR, comm);

        return recvbuf;
    }


    /**
     * Execute an all-gather communication with data of the given scalar type.
     * The returned vector contains the value provided by process j at int
     * j-th index.
     */
    template <typename T>
    std::vector<T> all_gather(const T& value) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");

        auto recvbuf = std::vector<T>(size(), T());
        MPI_Allgather(&value, sizeof(T), MPI_CHAR, &recvbuf[0], sizeof(T), MPI_CHAR, comm);
        return recvbuf;
    }


    /**
     * Execute an all-gather-v communication. This is a generalization of the
     * above, where each rank broadcasts to all others a container of items.
     * The size of the container to be broadcasted need not be the same on
     * every rank. The container broadcasted by rank j is returned in the j-th
     * index of the vector returned by this function.
     */
    template <typename T>
    std::vector<std::vector<T>> all_gather(const std::vector<T>& sendbuf) const
    {
        static_assert(std::is_trivially_copyable<T>::value, "type is not trivially copyable");

        auto recvcounts = all_gather(int(sendbuf.size() * sizeof(T)));
        auto recvdispls = std::vector<int>{0};
        std::partial_sum(recvcounts.begin(), recvcounts.end(), std::back_inserter(recvdispls));

        auto recvbuf = std::vector<T>(recvdispls.back() / sizeof(T));

        MPI_Allgatherv(
            &sendbuf[0], sendbuf.size() * sizeof(T), MPI_CHAR,
            &recvbuf[0], &recvcounts[0], &recvdispls[0], MPI_CHAR, comm);

        auto res = std::vector<std::vector<T>>(size());
        auto recv = recvbuf.begin();

        for (int i = 0; i < res.size(); ++i)
        {
            res[i].resize(recvcounts[i] / sizeof(T));

            for (int j = 0; j < res[i].size(); ++j)
            {
                res[i][j] = *recv++;
            }
        }
        return res;
    }


    // NOTE: The functions below must be generalized for other data types.
    // ========================================================================


    double reduce(std::size_t root, double local_value, operation op)
    {
        double result_value;
        MPI_Reduce(&local_value, &result_value, 1, MPI_DOUBLE, get_op(op), root, comm);
        return result_value;
    }

    double all_reduce(double local_value, operation op)
    {
        double result_value;
        MPI_Allreduce(&local_value, &result_value, 1, MPI_DOUBLE, get_op(op), comm);
        return result_value;
    }

    int all_reduce(int local_value, operation op)
    {
        int result_value;
        MPI_Allreduce(&local_value, &result_value, 1, MPI_INT, get_op(op), comm);
        return result_value;
    }


    //=========================================================================
    template<typename FunctionType, typename... Args>
    void invoke(FunctionType function, Args... args)
    {
        for (int i = 0; i < size(); ++i)
        {
            if (i == rank())
            {
                std::invoke(function, args...);
            }
            barrier();
        }
    }


private:
    // ========================================================================
    friend Communicator comm_world();
    MPI_Comm comm = MPI_COMM_NULL;
    bool owned = true;
};




// ============================================================================
mpi::Communicator mpi::comm_world()
{
    Communicator res;
    res.comm = MPI_COMM_WORLD;
    res.owned = false;
    return res;
}
