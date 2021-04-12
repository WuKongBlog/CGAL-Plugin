# CGAL-Plugin
## Fuzzy_rotate_cylinder
Find points inside a 3d rotate cylinder
Input:
- two endpoints of a line segment 
- radius of the cylinder
```C++
    // A simple example
    Point_d beginPt(58.5445, 90.3174, 149.704);
    Point_d endPt(58.1525, 92.8705, 149.728);
    std::vector<Point_d> vec_pt3;
    Traits::FT radius = 0.3;
    Fuzzy_rot_cylinder frc(beginPt, endPt, radius, 0.0);
    tree.search(std::back_inserter(vec_pt3), frc);
```
