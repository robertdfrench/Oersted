#ifndef OERSTED_DOUBLE2_H
#define OERSTED_DOUBLE2_H

struct double2 {
    double X;
    double Y;

    bool operator==(double2 const &d) const { return (X == d.X) && (Y == d.Y); };

    bool operator!=(double2 const &d) const { return (X != d.X) && (Y != d.Y); };

    bool operator>(double2 const &d) const { return (X > d.X) || (X == d.X && Y > d.Y); };

    bool operator<(double2 const &d) const { return (X < d.X) || (X == d.X && Y > d.Y); };
};

#endif // OERSTED_DOUBLE2_H
