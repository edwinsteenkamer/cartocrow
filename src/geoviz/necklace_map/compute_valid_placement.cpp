/*
The Necklace Map library implements the algorithmic geo-visualization
method by the same name, developed by Bettina Speckmann and Kevin Verbeek
at TU Eindhoven (DOI: 10.1109/TVCG.2010.180 & 10.1142/S021819591550003X).
Copyright (C) 2019  Netherlands eScience Center and TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Created by tvl (t.vanlankveld@esciencecenter.nl) on 07-01-2020
*/

#include "compute_valid_placement.h"

#include <glog/logging.h>

#include "geoviz/necklace_map/range.h"
#include "geoviz/necklace_map/detail/cycle_node.h"


namespace geoviz
{
namespace necklace_map
{
namespace detail
{

inline Number DistanceOnCircle(const Number& from_rad, const Number& to_rad)
{
  const Number dist = std::abs(to_rad - from_rad);
  return std::min(dist, M_2xPI - dist);
}

} // namespace detail

/**@class ComputeValidPlacement
 * @brief A functor to compute a valid placement for a collection of necklace beads.
 *
 * A placement for a set of necklace beads is a set of angles that describes the position of each bead on the necklace. A placement is valid if all scaled beads are inside their feasible interval and the distance between any two beads is at least some non-negative buffer distance. Note that this buffer distance guarantees that the beads do not overlap.
 *
 * There is often a range of valid placements. In this case, the placement is guided by an attraction-repulsion force: the beads are attracted to the center of their interval and repelled by the neighboring beads.
 *
 * Note that the placements are computed independently per necklace. This means that if a map contains multiple necklaces, no guarantees can be given about overlap of beads on different necklaces. However, such overlap can often be prevented manually by tuning the attraction-replusion force and the buffer distance.
 */

/**@brief Construct a new valid placement computation functor.
 * @param parameters the parameters describing the desired type of functor.
 * @return a unique pointer containing a new functor or a nullptr if the functor could not be constructed.
 */
ComputeValidPlacement::Ptr ComputeValidPlacement::New(const Parameters& parameters)
{
  switch (parameters.order_type)
  {
    case OrderType::kFixed:
      return Ptr(new ComputeValidPlacementFixedOrder(parameters.aversion_ratio, parameters.buffer_rad));
    default:
      return nullptr;
  }
}

ComputeValidPlacement::ComputeValidPlacement
(
  const Number& aversion_ratio,
  const Number& buffer_rad /*= 0*/
) :
  buffer_rad(buffer_rad),
  aversion_ratio(aversion_ratio)
{}

void ComputeValidPlacement::operator()(const Number& scale_factor, Necklace::Ptr& necklace) const
{
  const Number necklace_radius = necklace->shape->ComputeRadius();
  for (const Bead::Ptr& bead : necklace->beads)
  {
    // Compute the scaled covering radius.
    CHECK_NOTNULL(bead);
    const Number radius_scaled = scale_factor * bead->radius_base;
    bead->covering_radius_scaled_rad = std::asin(radius_scaled / necklace_radius);
  }

  // Compute the valid intervals.
  detail::ValidateScaleFactor validate(scale_factor, buffer_rad);
  const bool valid = validate(necklace);

  if (!valid)
    return;

  const double precision = 1e-7;

  const size_t num_beads = necklace->beads.size();
  for (int epoch = 0; epoch < 30; ++epoch)
  {
    for (size_t index_bead = 0; index_bead < num_beads; ++index_bead)
    {
      Bead::Ptr& bead = necklace->beads[index_bead];
      CHECK_NOTNULL(bead);

      const size_t index_prev = (index_bead + num_beads - 1) % num_beads;
      const size_t index_next = (index_bead + 1) % num_beads;

      Bead::Ptr& prev = necklace->beads[index_prev];
      const Bead::Ptr& next = necklace->beads[index_next];

      const Number offset_from_prev_rad = CircleRange(prev->angle_rad, bead->angle_rad).ComputeLength();
      const Number distance_neighbors_rad = CircleRange(prev->angle_rad, next->angle_rad).ComputeLength();
      const Number offset_from_centroid_rad = CircleRange
      (
        bead->feasible->ComputeCentroid(),
        bead->angle_rad
      ).ComputeLength();

      const Number& radius_bead = bead->covering_radius_scaled_rad;
      const Number& radius_prev = prev->covering_radius_scaled_rad;
      const Number& radius_next = next->covering_radius_scaled_rad;

      // Note that we cannot guarantee better than double precision because of the trigonometric functions.

      const double distance_to_prev_min = radius_prev + radius_bead + buffer_rad;
      const double distance_to_prev_max =
        (index_prev == index_next ? M_2xPI : distance_neighbors_rad) - radius_next - radius_bead - buffer_rad;

      // The 'bubble' is the largest range centered on the bead that does not contain the centroid.
      const double offset_prev_to_bubble =
        offset_from_centroid_rad < M_PI
        ? (offset_from_prev_rad - offset_from_centroid_rad)
        : (offset_from_prev_rad + (M_2xPI - offset_from_centroid_rad));

      const double w_0 =
        centroid_ratio * offset_prev_to_bubble * distance_to_prev_min * distance_to_prev_max -
        aversion_ratio * (distance_to_prev_min + distance_to_prev_max);
      const double w_1 =
        aversion_ratio * 2 -
        centroid_ratio *
        (
          (distance_to_prev_min + offset_prev_to_bubble) * (distance_to_prev_max + offset_prev_to_bubble) -
          offset_prev_to_bubble * offset_prev_to_bubble
        );
      const double w_2 = centroid_ratio * (distance_to_prev_min + distance_to_prev_max + offset_prev_to_bubble);
      const double w_3 = -centroid_ratio;

      // Solve w_3 * x^3 + w_2 * x^2 + w_1 * x + w_0 = 0 up to the specified precision.
      if (std::abs(w_3) < precision && std::abs(w_2) < precision)
      {
        const Number x = CircleRange::Modulo(-w_0 / w_1 + prev->angle_rad);

        if (!bead->feasible->IntersectsRay(x))
        {
          if (0 < 2 * offset_from_prev_rad - distance_to_prev_min + distance_to_prev_max)
            bead->angle_rad = bead->feasible->angle_cw_rad();
          else
            bead->angle_rad = bead->feasible->angle_ccw_rad();
        }
        else
          bead->angle_rad = x;
      }
      else
      {
        const double q = (3 * w_3 * w_1 - w_2 * w_2) / (9 * w_3 * w_3);
        const double r = (9 * w_3 * w_2 * w_1 - 27 * w_3 * w_3 * w_0 - 2 * w_2 * w_2 * w_2) / (54 * w_3 * w_3 * w_3);

        const double rho = std::max(std::sqrt(-q * q * q), std::abs(r));

        const double theta_3 = std::acos(r / rho) / 3;
        const double rho_3 = std::pow(rho, 1 / 3.0);

        const Number x = CircleRange::Modulo
        (
          prev->angle_rad -
          rho_3 * std::cos(theta_3) - w_2 / (3 * w_3) +
          rho_3 * std::sqrt(3.0) * std::sin(theta_3)
        );

        if (!bead->feasible->IntersectsRay(x))
        {
          if
          (
            0 <
            aversion_ratio * (2 * offset_from_prev_rad - (distance_to_prev_min + distance_to_prev_max)) +
            centroid_ratio * (offset_prev_to_bubble - offset_from_prev_rad) *
            (offset_from_prev_rad - distance_to_prev_min) * (offset_from_prev_rad - distance_to_prev_max)
          )
            bead->angle_rad = bead->feasible->angle_cw_rad();
          else
            bead->angle_rad = bead->feasible->angle_ccw_rad();
        }
        else
          bead->angle_rad = x;
      }
    }

    SwapBeads(necklace);
  }
}

void ComputeValidPlacement::operator()(const Number& scale_factor, std::vector<Necklace::Ptr>& necklaces) const
{
  for (Necklace::Ptr& necklace : necklaces)
    (*this)(scale_factor, necklace);
}


ComputeValidPlacementFixedOrder::ComputeValidPlacementFixedOrder
(
  const Number& aversion_ratio,
  const Number& buffer_rad /*= 0*/
) :
  ComputeValidPlacement(aversion_ratio, buffer_rad)
{}


ComputeValidPlacementAnyOrder::ComputeValidPlacementAnyOrder
(
  const Number& aversion_ratio,
  const Number& min_separation /*= 0*/
) :
  ComputeValidPlacement(aversion_ratio, min_separation) {}

void ComputeValidPlacementAnyOrder::SwapBeads(Necklace::Ptr& necklace) const
{
  const Number necklace_radius = necklace->shape->ComputeRadius();

  const size_t num_beads = necklace->beads.size();
  for (int index_bead = 0; index_bead < num_beads; index_bead++)
  {
    const size_t index_next = (index_bead+1)%num_beads;

    Bead::Ptr& bead = necklace->beads[index_bead];
    Bead::Ptr& next = necklace->beads[index_next];

    const Number& radius_bead = bead->covering_radius_scaled_rad;
    const Number& radius_next = next->covering_radius_scaled_rad;

    // Note that for the swapped angles, the buffers cancel each other out.
    const Number swapped_angle_bead_rad = CircleRange::Modulo(next->angle_rad + radius_next - radius_bead);
    const Number swapped_angle_next_rad = CircleRange::Modulo(bead->angle_rad - radius_bead + radius_next);

    if (bead->feasible->IntersectsRay(swapped_angle_bead_rad) && next->feasible->IntersectsRay(swapped_angle_next_rad))
    {
      const Number centroid_bead_rad = bead->feasible->ComputeCentroid();
      const Number centroid_next_rad = next->feasible->ComputeCentroid();

      const Number dist_original_bead = detail::DistanceOnCircle(bead->angle_rad, centroid_bead_rad);
      const Number dist_original_next = detail::DistanceOnCircle(next->angle_rad, centroid_next_rad);
      const Number dist_swapped_bead = detail::DistanceOnCircle(swapped_angle_bead_rad, centroid_bead_rad);
      const Number dist_swapped_next = detail::DistanceOnCircle(swapped_angle_next_rad, centroid_next_rad);

      const Number cost_original = dist_original_bead * dist_original_bead + dist_original_next * dist_original_next;
      const Number cost_swapped = dist_swapped_bead * dist_swapped_bead + dist_swapped_next * dist_swapped_next;

      if (cost_swapped < cost_original)
      {
        // Swap the beads.
        bead->angle_rad = swapped_angle_bead_rad;
        next->angle_rad = swapped_angle_next_rad;

        std::swap(necklace->beads[index_bead], necklace->beads[index_next]);
      }
    }
  }
}

} // namespace necklace_map
} // namespace geoviz