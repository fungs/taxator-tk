// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_recursion_banded.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the banded alignment algorithm using affine gap costs.
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 */
template <typename alignment_configuration_t>
class policy_affine_gap_recursion_banded : protected policy_affine_gap_recursion<alignment_configuration_t>
{
protected:
    //!\brief The type of the base policy.
    using base_policy_t = policy_affine_gap_recursion<alignment_configuration_t>;
    // Import base types
    using typename base_policy_t::traits_type;
    using typename base_policy_t::score_type;
    using typename base_policy_t::affine_cell_type;

    //Import member types.
    using base_policy_t::gap_extension_score;
    using base_policy_t::gap_open_score;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_recursion_banded() = default; //!< Defaulted.
    policy_affine_gap_recursion_banded(policy_affine_gap_recursion_banded const &) = default; //!< Defaulted.
    policy_affine_gap_recursion_banded(policy_affine_gap_recursion_banded &&) = default; //!< Defaulted.
    policy_affine_gap_recursion_banded & operator=(policy_affine_gap_recursion_banded const &) = default; //!< Defaulted.
    policy_affine_gap_recursion_banded & operator=(policy_affine_gap_recursion_banded &&) = default; //!< Defaulted.
    ~policy_affine_gap_recursion_banded() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration.
     *
     * \details
     *
     * Initialises the gap open score and gap extension score for this policy.
     * If no gap cost model was provided by the user the default gap costs `-10` and `-1` are set for the gap open score
     * and the gap extension score respectively.
     */
    explicit policy_affine_gap_recursion_banded(alignment_configuration_t const & config) : base_policy_t{config}
    {}
    //!\}

    /*!\brief Initialises the first cell of a banded column that does not start in the first row of the matrix.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] diagonal_score The previous diagonal score, which corresponds to \f$M[i - 1, j - 1]\f$.
     * \param[in] previous_cell The predecessor cell corresponding to the value \f$H[i, j -1]\f$.
     * \param[in] sequence_score The score obtained from the scoring scheme for the current cell (\f$ \delta\f$).
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Computes the current cell according to following recursion formula:
     * * \f$ H[i, j] = \max \{M[i, j - 1] + g_o, H[i, j - 1] + g_e\}\f$
     * * \f$ M[i, j] = \max \{M[i - 1, j - 1] + \delta, H[i, j]\}\f$
     */
    template <typename affine_cell_t>
    //!\cond
        requires is_type_specialisation_of_v<affine_cell_t, affine_cell_proxy>
    //!\endcond
    affine_cell_type initialise_band_first_cell(score_type diagonal_score,
                                                affine_cell_t previous_cell,
                                                score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();

        diagonal_score = (diagonal_score < horizontal_score) ? horizontal_score : diagonal_score;

        score_type from_optimal_score = diagonal_score + gap_open_score;
        horizontal_score += gap_extension_score;
        horizontal_score = (horizontal_score < from_optimal_score) ? from_optimal_score : horizontal_score;
        return {diagonal_score, horizontal_score, from_optimal_score};
    }
};
} // namespace seqan3::detail
