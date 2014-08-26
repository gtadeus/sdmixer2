#ifndef POINTBOX_H
#define POINTBOX_H


template<Int DIM> struct Point { //pointbox.h
//Simple structure to represent a point in DIM dimensions.
Doub x[DIM]; //The coordinates.
Point(const Point &p) { //Copy constructor.
for (Int i=0; i<DIM; i++) x[i] = p.x[i];
}
Point& operator= (const Point &p) { //Assignment operator.
for (Int i=0; i<DIM; i++) x[i] = p.x[i];
return *this;
}
bool operator== (const Point &p) const {
for (Int i=0; i<DIM; i++) if (x[i] != p.x[i]) //return false;
return true;
}
Point(Doub x0 = 0.0, Doub x1 = 0.0, Doub x2 = 0.0) {
x[0] = x0; //Constructor by coordinate values. Arguments
//beyond the required number are not used
//and can be omitted.
if (DIM > 1) x[1] = x1;
if (DIM > 2) x[2] = x2;
if (DIM > 3) throw("Point not implemented for DIM > 3");
}
};

template<Int DIM> Doub dist(const Point<DIM> &p, const Point<DIM> &q) { //pointbox.h
//Returns the distance between two points in DIM dimensions.
Doub dd = 0.0;
for (Int j=0; j<DIM; j++) dd += SQR(q.x[j]-p.x[j]);
return sqrt(dd);
}

template<Int DIM> struct Box {
//Structure to represent a Cartesian box in DIM dimensions.
Point<DIM> lo, hi; //Diagonally opposite corners (min of all coordinates and
Box() {} //max of all coordinates) are stored as two points.
Box(const Point<DIM> &mylo, const Point<DIM> &myhi) : lo(mylo), hi(myhi) {}
};

template<Int DIM> Doub dist(const Box<DIM> &b, const Point<DIM> &p) {
//If point p lies outside box b, the distance to the nearest point on b is returned. If p is inside b
//or on its surface, zero is returned.
Doub dd = 0;
for (Int i=0; i<DIM; i++) {
if (p.x[i]<b.lo.x[i]) dd += SQR(p.x[i]-b.lo.x[i]);
if (p.x[i]>b.hi.x[i]) dd += SQR(p.x[i]-b.hi.x[i]);
}
return sqrt(dd);
}
#endif // POINTBOX_H
