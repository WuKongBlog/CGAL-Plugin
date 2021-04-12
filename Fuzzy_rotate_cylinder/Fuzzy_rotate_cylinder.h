// Edit by Wsh 2021-04-12
// Motivated by
// https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
// https://www.gamedev.net/forums/topic/338522-bounding-box-for-a-cylinder/
// 

#ifndef CGAL_FUZZY_ROTATE_CYLINDER_H
#define CGAL_FUZZY_ROTATE_CYLINDER_H

#include <CGAL/license/Spatial_searching.h>

#include <algorithm>

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Search_traits_adapter.h>

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/enable_if.hpp>

namespace CGAL
{

	namespace internal
	{
		template <class SearchTraits, class Point>
		struct wshIs_from_point_from_adapter_traits
		{
			typedef boost::false_type type;
		};


		template <class K, class PM, class Base, class Point>
		struct wshIs_from_point_from_adapter_traits<Search_traits_adapter<K, PM, Base>, Point>
		{
			typedef typename boost::is_same<Point, typename Base::Point_d> type;
		};
	} //namespace internal

	template <class SearchTraits>
	class Fuzzy_rot_cylinder
	{
		SearchTraits traits;

	public:

		typedef typename SearchTraits::Point_d Point_d;
		typedef typename SearchTraits::Iso_box_d Iso_box_d;
		typedef typename SearchTraits::FT FT;
		typedef typename SearchTraits::Dimension Dimension;
		typedef typename SearchTraits::Construct_min_vertex_d Construct_min_vertex_d;
		typedef typename SearchTraits::Construct_max_vertex_d Construct_max_vertex_d;
		typedef typename SearchTraits::Cartesian_const_iterator_d Cartesian_const_iterator_d;
		typedef typename SearchTraits::Construct_cartesian_const_iterator_d Construct_cartesian_const_iterator_d;


	private:

		typename boost::remove_cv<
			typename boost::remove_reference< typename Construct_min_vertex_d::result_type >::type
		>::type min, max;
		Cartesian_const_iterator_d min_begin, max_begin, cyl_p_begin, cyl_q_begin;
		FT eps, rad;
		unsigned int dim;

		//constructor implementation
		template <class Point, class Construct_iso_box_d>
		void construct(const Point& p, const Point& q, FT radius)
		{
			Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
			Cartesian_const_iterator_d begin = construct_it(p),
				end = construct_it(p, 1);
			dim = static_cast<unsigned int>(end - begin);

			// only support 3d
			CGAL_precondition(dim == 3);

			Point maxpt, minpt;

			Cartesian_const_iterator_d
				pit = construct_it(p),
				qit = construct_it(q);

			//bounding box for a cylinder
			for (unsigned int i = 0; i < dim; ++pit, ++qit, ++i)
			{
				if ((*pit) < (*qit))
				{
					minpt.at(i) = (*pit) - radius;
					maxpt.at(i) = (*qit) + radius;
				}
				else
				{
					minpt.at(i) = (*qit) - radius;
					maxpt.at(i) = (*pit) + radius;
				}
			}

			Iso_box_d box = Construct_iso_box_d()(minpt, maxpt);
			Construct_min_vertex_d construct_min_vertex_d;
			Construct_max_vertex_d construct_max_vertex_d;
			min = construct_min_vertex_d(box);
			max = construct_max_vertex_d(box);
			min_begin = construct_it(min);
			max_begin = construct_it(max);

			cyl_p_begin = construct_it(p);
			cyl_q_begin = construct_it(q);
		}

	public:

		// default constructor
		Fuzzy_rot_cylinder(const SearchTraits& traits_ = SearchTraits()) :traits(traits_) {}

		// constructor
		Fuzzy_rot_cylinder(const Point_d& p, const Point_d& q, FT radius_, FT epsilon_ = FT(0), const SearchTraits& traits_ = SearchTraits())
			: traits(traits_), eps(epsilon_), rad(radius_)
		{
			CGAL_precondition(epsilon_ == 0.0);
			construct<Point_d, typename SearchTraits::Construct_iso_box_d>(p, q, radius_);
		}

		//additional constructor if SearchTraits = Search_traits_adapter
		template <class Point>
		Fuzzy_rot_cylinder(const Point& p, const Point& q, FT epsilon = FT(0), const SearchTraits& traits_ = SearchTraits(),
			typename boost::enable_if<typename internal::wshIs_from_point_from_adapter_traits<SearchTraits, Point>::type>::type* = 0)
			: traits(traits_), eps(epsilon)
		{
			CGAL_precondition(epsilon == 0.0);
			construct<Point, typename SearchTraits::Construct_iso_box_d>(p, q);
		}

		bool contains(const Point_d& p) const
		{
			Construct_cartesian_const_iterator_d construct_it = traits.construct_cartesian_const_iterator_d_object();
			Cartesian_const_iterator_d pit = construct_it(p);
			Cartesian_const_iterator_d cyl_pit = cyl_p_begin, cyl_qit = cyl_q_begin;
			
			//check if a 3d point is inside a cylinder
			FT
				pt0 = *pit,
				pt1 = *(pit + 1),
				pt2 = *(pit + 2);

			FT
				p0 = *cyl_pit,
				p1 = *(cyl_pit + 1),
				p2 = *(cyl_pit + 2);

			FT
				q0 = *cyl_qit,
				q1 = *(cyl_qit + 1),
				q2 = *(cyl_qit + 2);

			FT 
				check1 = (pt0 - p0) * (q0 - p0) + (pt1 - p1) * (q1 - p1) + (pt2 - p2) * (q2 - p2),
				check2 = (pt0 - q0) * (q0 - p0) + (pt1 - q1) * (q1 - p1) + (pt2 - q2) * (q2 - p2);

			if (check1 < 0 || check2 > 0)
			{
				return false;
			}

			FT
				temp0 = (pt1 - p1) * (q2 - p2) - (pt2 - p2) * (q1 - p1),
				temp1 = (pt2 - p2) * (q0 - p0) - (pt0 - p0) * (q2 - p2),
				temp2 = (pt0 - p0) * (q1 - p1) - (pt1 - p1) * (q0 - p0);

			FT
				check3 = CGAL::sqrt(temp0 * temp0 + temp1 * temp1 + temp2 * temp2),
				check4 = rad * CGAL::sqrt((q0 - p0) * (q0 - p0) + (q1 - p1) * (q1 - p1) + (q2 - p2) * (q2 - p2));

			if (check3 > check4)
			{
				return false;
			}

			return true;
		}

		bool inner_range_intersects(const Kd_tree_rectangle<FT, Dimension>& rectangle) const
		{
			// test whether the box eroded by 'eps' intersects 'rectangle'
			Cartesian_const_iterator_d minit = min_begin, maxit = max_begin;
			for (unsigned int i = 0; i < dim; ++i, ++minit, ++maxit) 
			{
				if (((*maxit) - eps < rectangle.min_coord(i))
					|| ((*minit) + eps > rectangle.max_coord(i)))
					return false;
			}
			return true;
		}

		bool outer_range_contains(const Kd_tree_rectangle<FT, Dimension>& rectangle) const
		{
			// Cause It's hard to determine if a iso-box in a rot-cylinder and usually the radius is small
			// So here just make sure every point will do a check
			return false;
		}
	}; // class Fuzzy_rot_cylinder

} // namespace CGAL
#endif // Fuzzy_rot_cylinder_H
