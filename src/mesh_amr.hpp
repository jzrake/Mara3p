/**
 ==============================================================================
 Copyright 2020, Jonathan Zrake

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
#include "core_ndarray.hpp"
#include "core_ndarray_ops.hpp"
#include "core_bsp_tree.hpp"
#include "core_bqo_tree.hpp"
#include "core_numeric_array.hpp"
#include "parallel_computable.hpp"




//=============================================================================
namespace amr
{




//=============================================================================
struct tile_blocks
{
    template<typename P>
    auto operator()(nd::array_t<P, 1> blocks) const
    {
        static_assert(nd::is_array<std::invoke_result_t<P, nd::uivec_t<1>>>::value,
            "must be a 1d array of 1d arrays");

        for (auto block : blocks)
            if (shape(block) != shape(front(blocks)))
                throw std::invalid_argument("amr::tile_blocks (supplied blocks have non-uniform shapes)");

        auto bi = shape(blocks, 0);
        auto ni = shape(front(blocks), 0);

        return nd::make_array(nd::indexing([ni, blocks] (auto i)
        {
            auto b = nd::uivec(i / ni);
            auto c = nd::uivec(i % ni);
            return blocks(b)(c);
        }), nd::uivec(bi * ni));
    }

    template<typename P>
    auto operator()(nd::array_t<P, 2> blocks) const
    {
        static_assert(nd::is_array<std::invoke_result_t<P, nd::uivec_t<2>>>::value,
            "must be a 2d array of 2d arrays");

        for (auto block : blocks)
            if (shape(block) != shape(front(blocks)))
                throw std::invalid_argument("amr::tile_blocks (supplied blocks have non-uniform shapes)");

        auto bi = shape(blocks, 0);
        auto bj = shape(blocks, 1);
        auto ni = shape(front(blocks), 0);
        auto nj = shape(front(blocks), 1);

        return nd::make_array(nd::indexing([ni, nj, blocks] (auto i, auto j)
        {
            auto b = nd::uivec(i / ni, j / nj);
            auto c = nd::uivec(i % ni, j % nj);
            return blocks(b)(c);
        }), nd::uivec(bi * ni, bj * nj));
    }
};




//=============================================================================
struct downsample
{
    template<typename P>
    auto operator()(nd::array_t<P, 1> a) const
    {
        auto [ni] = as_tuple(shape(a));

        if (ni % 2 != 0)
        {
            throw std::invalid_argument("amr::downsample (array size must be even on all axes)");
        }

        return nd::make_array(nd::indexing([a] (auto i)
        {
            return 0.5 * (
                a(i * 2 + 0) +
                a(i * 2 + 1));
        }), nd::uivec(ni / 2));
    }

    template<typename P>
    auto operator()(nd::array_t<P, 2> a) const
    {
        auto [ni, nj] = as_tuple(shape(a));

        if (ni % 2 != 0 || nj % 2 != 0)
        {
            throw std::invalid_argument("amr::downsample (array size must be even on all axes)");
        }

        return nd::make_array(nd::indexing([a] (auto i, auto j)
        {
            return 0.25 * (
                a(i * 2 + 0, j * 2 + 0) +
                a(i * 2 + 0, j * 2 + 1) +
                a(i * 2 + 1, j * 2 + 0) +
                a(i * 2 + 1, j * 2 + 1));
        }), nd::uivec(ni / 2, nj / 2));
    }
};




//=============================================================================
struct upsample
{
    template<typename P>
    auto operator()(nd::array_t<P, 1> a)
    {
        auto [ni] = as_tuple(shape(a));

        return nd::make_array(nd::indexing([a] (auto i)
        {
            return a(i / 2);
        }), nd::uivec(ni * 2));
    }

    template<typename P>
    auto operator()(nd::array_t<P, 2> a)
    {
        auto [ni, nj] = as_tuple(shape(a));

        return nd::make_array(nd::indexing([a] (auto i, auto j)
        {
            return a(i / 2, j / 2);
        }), nd::uivec(ni * 2, nj * 2));
    }
};




//=============================================================================
struct coarsen_blocks
{
    template<typename P>
    auto operator()(nd::array_t<P, 1> v0, nd::array_t<P, 1> v1) const
    {
        return nd::make_array(nd::indexing([=] (auto i)
        {
            if (i == 0) return v0;
            if (i == 1) return v1;
            throw std::invalid_argument("amr::coarsen_blocks");
        }), nd::uivec(2)) | tile_blocks() | downsample() | nd::to_shared();
    }

    template<typename P>
    auto operator()(nd::array_t<P, 2> v00, nd::array_t<P, 2> v10, nd::array_t<P, 2> v01, nd::array_t<P, 2> v11) const
    {
        return nd::make_array(nd::indexing([=] (auto i, auto j)
        {
            if (i == 0 && j == 0) return v00;
            if (i == 0 && j == 1) return v01;
            if (i == 1 && j == 0) return v10;
            if (i == 1 && j == 1) return v11;
            throw std::invalid_argument("amr::coarsen_blocks");
        }), nd::uivec(2, 2)) | tile_blocks() | downsample() | nd::to_shared();
    }
};




//=============================================================================
struct split_block
{
    split_block(unsigned which) : which(which) {}

    template<typename P>
    auto operator()(nd::array_t<P, 1> block) const
    {
        auto A = [block] (nd::uint a)
        {
            return nd::indexing([block, a] (nd::uint i)
            {
                return block(a + i);
            });
        };

        auto [ni] = as_tuple(shape(block));

        switch (which)
        {
            case 0: return nd::make_array(A(     0), nd::uivec(ni / 2));
            case 1: return nd::make_array(A(ni / 2), nd::uivec(ni / 2));
        }
        throw std::invalid_argument("amr::split_block");
    }

    template<typename P>
    auto operator()(nd::array_t<P, 2> block) const
    {
        auto A = [block] (nd::uint a, nd::uint b)
        {
            return nd::indexing([block, a, b] (nd::uint i, nd::uint j)
            {
                return block(a + i, b + j);
            });
        };

        auto [ni, nj] = as_tuple(shape(block));

        switch (which)
        {
            case 0: return nd::make_array(A(     0,      0), nd::uivec(ni / 2, nj / 2));
            case 1: return nd::make_array(A(ni / 2,      0), nd::uivec(ni / 2, nj / 2));
            case 2: return nd::make_array(A(     0, nj / 2), nd::uivec(ni / 2, nj / 2));
            case 3: return nd::make_array(A(ni / 2, nj / 2), nd::uivec(ni / 2, nj / 2));
        }
        throw std::invalid_argument("amr::split_block");
    }
    unsigned which = 0;
};




//=============================================================================
struct refine_block
{
    refine_block(unsigned which) : which(which) {}

    template<typename ValueType, unsigned long Rank>
    auto operator()(nd::shared_array<ValueType, Rank> v) const
    {
        return v | upsample() | split_block(which) | nd::to_shared();
    }
    unsigned which = 0;
};




//=============================================================================
template<typename ValueType>
auto collapse(bsp::shared_tree<mpr::computable<nd::shared_array<ValueType, 1>>, 2> tree)
{
    return bsp::collapse(tree, [] (auto blocks)
    {
        return mpr::zip(blocks[0], blocks[1]) | mpr::mapv(coarsen_blocks());
    });
}




//=============================================================================
template<typename ValueType>
auto collapse(bsp::shared_tree<mpr::computable<nd::shared_array<ValueType, 2>>, 4> tree)
{
    return bsp::collapse(tree, [] (auto blocks)
    {
        return mpr::zip(blocks[0], blocks[1], blocks[2], blocks[3]) | mpr::mapv(coarsen_blocks());
    });
}




//=============================================================================
template<typename ArrayType, unsigned long Ratio>
auto get_or_create_block(bsp::shared_tree<mpr::computable<ArrayType>, Ratio> tree, mesh::block_index_t<bsp::log2<Ratio>::value> block)
{
    // If the tree has a value at the target block, then return that value.

    if (contains(tree, block))
    {
        return value_at(tree, block) | mpr::map([] (auto v) { return v; });
    }

    // If the tree has a value at the node above the target block, then
    // refine the data on that node and select the array in the block's
    // orthant.

    if (contains(tree, parent_index(block)))
    {
        return value_at(tree, parent_index(block)) | mpr::map(refine_block(to_integral(orthant(relative_to_parent(block)))));
    }

    // If the target block is not a leaf, then tile and downsample its child
    // blocks.

    if (contains_node(tree, block))
    {
        return collapse(node_at(tree, block));
    }

    throw std::logic_error("amr::get_or_create_block (invalid mesh topology)");
};




//=============================================================================
template<typename AtFunction>
auto extend_block(AtFunction at, mesh::block_index_t<1> block)
{
    auto c0 = at(prev_on(block, 0));
    auto c1 = at(block);
    auto c2 = at(next_on(block, 0));

    return mpr::zip(c0, c1, c2)
    | mpr::mapv([] (auto c0, auto c1, auto c2)
    {
        auto nx = shape(c1, 0);
        auto cs = std::array{c0, c1, c2};

        return nd::make_array(nd::indexing([cs, nx] (auto i)
        {
            auto bi = i < 2 ? 0 : (i >= nx + 2 ? 2 : 1);
            auto ii = i < 2 ? i - 2 + nx : (i >= nx + 2 ? i - 2 - nx : i - 2);
            return cs[bi](ii);
        }), nd::uivec(nx + 4)) | nd::to_shared();
    });
}




//=============================================================================
template<typename AtFunction>
auto extend_block(AtFunction at, mesh::block_index_t<2> block)
{
    auto c11 = at(block);
    auto c00 = at(prev_on(prev_on(block, 0), 1));
    auto c02 = at(prev_on(next_on(block, 0), 1));
    auto c20 = at(next_on(prev_on(block, 0), 1));
    auto c22 = at(next_on(next_on(block, 0), 1));
    auto c01 = at(prev_on(block, 0));
    auto c21 = at(next_on(block, 0));
    auto c10 = at(prev_on(block, 1));
    auto c12 = at(next_on(block, 1));

    return mpr::zip(c00, c01, c02, c10, c11, c12, c20, c21, c22)
    | mpr::mapv([] (auto c00, auto c01, auto c02, auto c10, auto c11, auto c12, auto c20, auto c21, auto c22)
    {
        auto nx = shape(c11, 0);
        auto ny = shape(c11, 1);
        auto cs = std::array{
            std::array{c00, c01, c02},
            std::array{c10, c11, c12},
            std::array{c20, c21, c22},
        };

        return nd::make_array(nd::indexing([cs, nx, ny] (auto i, auto j)
        {
            auto bi = i < 2 ? 0 : (i >= nx + 2 ? 2 : 1);
            auto bj = j < 2 ? 0 : (j >= ny + 2 ? 2 : 1);
            auto ii = i < 2 ? i - 2 + nx : (i >= nx + 2 ? i - 2 - nx : i - 2);
            auto jj = j < 2 ? j - 2 + ny : (j >= ny + 2 ? j - 2 - ny : j - 2);
            return cs[bi][bj](ii, jj);
        }), nd::uivec(nx + 4, ny + 4)) | nd::to_shared();
    });
}




/**
 * @brief      Determine whether each an index in a tree has any over-refined
 *             neighbors. A node adjacent to a leaf node L is over-refined if
 *             its maximum depth is more than one greater than the depth of L.
 *
 * @param[in]  tree  The tree to check
 *
 * @return     True or false
 */
inline auto has_over_refined_neighbors(bsp::shared_tree<mesh::block_index_t<1>, 2> tree)
{
    return [tree] (mesh::block_index_t<1> block)
    {
        return
        (contains_node(tree, next_on(block, 0)) && depth(node_at(tree, next_on(block, 0))) > 1) ||
        (contains_node(tree, prev_on(block, 0)) && depth(node_at(tree, prev_on(block, 0))) > 1);
    };
}

inline auto has_over_refined_neighbors(bsp::shared_tree<mesh::block_index_t<2>, 4> tree)
{
    return [tree] (mesh::block_index_t<2> block)
    {
        return
        (contains_node(tree, next_on(block, 0)) && depth(node_at(tree, next_on(block, 0))) > 1) ||
        (contains_node(tree, prev_on(block, 0)) && depth(node_at(tree, prev_on(block, 0))) > 1) ||
        (contains_node(tree, next_on(block, 1)) && depth(node_at(tree, next_on(block, 1))) > 1) ||
        (contains_node(tree, prev_on(block, 1)) && depth(node_at(tree, prev_on(block, 1))) > 1);
    };
}




/**
 * @brief      Return a mesh tree (binary, quad, or oct) guaranteed not to have
 *             any over-refined neighbors. This is accomplished only by
 *             repeatedly branching nodes that have over-refined neighbors.
 *
 * @param[in]  tree  The tree of vertex blocks.
 *
 * @return     A tree without any over-refined neighbors
 */
template<unsigned long Ratio>
inline auto valid_mesh_tree(bsp::shared_tree<mesh::block_index_t<bsp::log2<Ratio>::value>, Ratio> tree)
{
    while (reduce(bsp::map(indexes(tree), has_over_refined_neighbors(tree)), std::logical_or<>(), false))
    {
        tree = bsp::branch_if(tree, mesh::child_indexes<bsp::log2<Ratio>::value>, has_over_refined_neighbors(tree));
    }
    return tree;
}




//=============================================================================
template<typename ValueType, unsigned long Ratio>
auto weighted_sum_tree(
    bsp::shared_tree<mpr::computable<ValueType>, Ratio> s,
    bsp::shared_tree<mpr::computable<ValueType>, Ratio> t, double b)
{
    return bsp::zip(s, t) | bsp::mapvs([b] (auto s, auto t)
    {
        return mpr::zip(s, t) | mpr::mapv([b] (auto s, auto t)
        {
            return nd::to_shared(s * b + t * (1.0 - b));
        }, "RK");
    });
}

} // namespace amr
