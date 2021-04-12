# CGAL-Plugin
woh
```C++
    // Example
    Point_d beginPt(58.5445, 90.3174, 149.704);
    Point_d endPt(58.1525, 92.8705, 149.728);
    std::vector<Point_d> vec_pt3;
    Traits::FT radius = 0.3;
    Fuzzy_rot_cylinder frc(beginPt, endPt, radius, 0.0);
    tree.search(std::back_inserter(vec_pt3), frc);
```
