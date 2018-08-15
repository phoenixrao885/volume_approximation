// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef BALLINTERSECTCONVEX_H
#define BALLINTERSECTCONVEX_H

// ball type
struct Ball{
public:
    Ball(Point c, NT R) : _c(c),	 _R(R) {}

    Point center(){
        return _c;
    }

    NT squared_radius(){
        return _R;
    }

    NT radius(){
        return std::sqrt(_R);
    }

    int is_in(Point p){
        if ((p-_c).squared_length() <= _R)
            return -1;
        else return 0;
    }

    std::pair<NT,NT> line_intersect(Point r,
                                          Point v){
        //Point::Cartesian_const_iterator rit;
        //rit=r.cartesian_begin();
        typename std::vector<NT>::iterator rit=r.iter_begin();
        //Point::Cartesian_const_iterator vit;
        //vit=v.cartesian_begin();
        typename std::vector<NT>::iterator vit=v.iter_begin();
        //Point::Cartesian_const_iterator cit;
        //cit=_c.cartesian_begin();
        typename std::vector<NT>::iterator cit=_c.iter_begin();
        Point rc = r - _c;
        //Vector::Cartesian_const_iterator rcit;
        //rcit=rc.cartesian_begin();
        typename std::vector<NT>::iterator rcit=rc.iter_begin();
        NT vrc(0);
        NT v2(0);
        NT rc2(0);
        for( ; cit < _c.iter_end() ; ++rcit, ++cit, ++rit, ++vit){
            vrc += *vit * (*rcit);
            v2 += *vit * (*vit);
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - _R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);
        //return std::pair<Point,Point> ((lamda1*v)+r,(lamda2*v)+r);
        return std::pair<NT,NT> (lamda1,lamda2);
    }

    std::pair<NT,NT> line_intersect_coord(Point r,
                                          int rand_coord){

        Point rc = r - _c;
        //Vector::Cartesian_const_iterator rcit;
        //rcit = rc.cartesian_begin();
        typename std::vector<NT>::iterator rcit=rc.iter_begin();
        NT vrc = *(rcit + rand_coord);

        NT v2 = NT(1);
        NT rc2(0);
        for( ; rcit < rc.iter_end() ; ++rcit){
            rc2 += *rcit * (*rcit);
        }

        NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - _R));
        NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);
        NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);

        return std::pair<NT,NT> (lamda1,lamda2);

    }

private:
    Point  _c; //center
    NT     _R; //SQUARED radius !!!
};


template <class T, typename FT>
class BallIntersectPolytope {
public:
    BallIntersectPolytope(T &P, Ball &B) : _P(P), _B(B) {};
    
    T first() { return _P; }
    Ball second() { return _B; }
    
    int is_in(Point p) {
        if (_B.is_in(p) == -1)
            return _P.is_in(p);
        return 0;
    }

    int num_of_hyperplanes(){
        return _P.num_of_hyperplanes();
    }

    int dimension(){
        return _P.dimension();
    }

    std::pair<FT,FT> line_intersect(Point r,
                                          Point v) {

        std::pair <FT, FT> polypair = _P.line_intersect(r, v);
        std::pair <FT, FT> ballpair = _B.line_intersect(r, v);
        return std::pair<FT, FT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    //First coordinate ray shooting intersecting convex body
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          int rand_coord,
                                          std::vector<FT> &lamdas) {

        std::pair <FT, FT> polypair = _P.line_intersect_coord(r, rand_coord, lamdas);
        std::pair <FT, FT> ballpair = _B.line_intersect_coord(r, rand_coord);
        return std::pair<FT, FT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    //Not the first coordinate ray shooting intersecting convex body
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<FT> &lamdas) {

        std::pair <FT, FT> polypair = _P.line_intersect_coord(r, r_prev, rand_coord, rand_coord_prev, lamdas);
        std::pair <FT, FT> ballpair = _B.line_intersect_coord(r, rand_coord);
        return std::pair<FT, FT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }

    std::pair<FT,FT> query_dual(Point &p, int rand_coord) {
        std::pair <FT, FT> polypair = _P.query_dual(p, rand_coord);
        std::pair <FT, FT> ballpair = _B.line_intersect_coord(p, rand_coord);
        return std::pair<FT, FT>(std::min(polypair.first, ballpair.first),
                                 std::max(polypair.second, ballpair.second));
    }
    
private:
    T    _P;
    Ball _B;
};


// Convex body defined as an intersection of balls
template <typename FT>
class IntersectionOfBalls {
private:
    std::vector<Ball> balls;
    int _d;
public:
    IntersectionOfBalls (int dim, std::vector<Ball> vecballs) : _d(dim), balls(vecballs) {};

    int dimension() {
        return _d;
    }

    int num_of_hyperplanes() {
        return 0;
    }

    void add_ball (Ball B) {
        balls.push_back(B);
    }

    int is_in(Point p) {
        typename std::vector<Ball>::iterator itB = balls.begin();
        for ( ; itB!=balls.end(); itB++) {
            if ((*itB).is_in(p)==0) {
                return 0;
            }
        }
        return -1;
    }

    std::pair<FT,FT> line_intersect(Point r,
                                    Point v) {
        FT min_plus = FT(maxNT), max_minus = FT(minNT);
        std::pair <FT, FT> ballpair;
        typename std::vector<Ball>::iterator itB = balls.begin();
        for ( ; itB!=balls.end(); itB++) {
            std::pair <FT, FT> ballpair = (*itB).line_intersect(r, v);
            if (ballpair.first < min_plus) min_plus = ballpair.first;
            if (ballpair.second > max_minus) max_minus = ballpair.first;
        }
        return std::pair<FT, FT> (min_plus, max_minus);
    }

    //First coordinate ray shooting intersecting convex body
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          int rand_coord,
                                          std::vector<FT> &lamdas) {
        FT min_plus = FT(maxNT), max_minus = FT(minNT);
        std::pair <FT, FT> ballpair;
        typename std::vector<Ball>::iterator itB = balls.begin();
        for ( ; itB!=balls.end(); itB++) {
            std::pair <FT, FT> ballpair = (*itB).line_intersect_coord(r, rand_coord);
            if (ballpair.first < min_plus) min_plus = ballpair.first;
            if (ballpair.second > max_minus) max_minus = ballpair.first;
        }
        return std::pair<FT, FT> (min_plus, max_minus);
    }

    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<FT> &lamdas) {
        FT min_plus = FT(maxNT), max_minus = FT(minNT);
        std::pair <FT, FT> ballpair;
        typename std::vector<Ball>::iterator itB = balls.begin();
        for ( ; itB!=balls.end(); itB++) {
            std::pair <FT, FT> ballpair = (*itB).line_intersect_coord(r, rand_coord);
            if (ballpair.first < min_plus) min_plus = ballpair.first;
            if (ballpair.second > max_minus) max_minus = ballpair.first;
        }
        return std::pair<FT, FT> (min_plus, max_minus);
    }


};


template <class T1 , class T2>
class PolytopeIntersectEllipsoid {
private:
    T1 P;
    T2 E;
    typedef typename T2::K 	K;
public:
    PolytopeIntersectEllipsoid(T1 &Pin, T2 &Ein) : P(Pin), E(Ein) {};
    
    T1 first() { return P; }
    T2 second() { return E; }
    
    int is_in(Point p){
        //std::cout << "calling is in"<<std::endl;
        if(P.is_in(p)==-1)
            return E.is_in(p);
        return 0;
    }
    
    int num_of_hyperplanes(){
        return P.num_of_hyperplanes();
    }

    int dimension(){
        return P.dimension();
    }
    
    std::pair<Point,Point> line_intersect(Point r,
                                          Point v){

        std::pair<Point,Point> polypair = P.line_intersect(r,v);
        std::pair<Point,Point> returnpair;
        std::pair<Point,Point> ellpair;
        bool ellinter=false;

        //check the first intersection point if it is inside ball
        if(E.is_in(polypair.first)){
            returnpair.first = polypair.first;
        }else{
            ellinter=true;
            //compute the intersection with ball
            ellpair = E.line_intersect(r,v);
            returnpair.first = ellpair.first;
        }
        //check the second intersection point
        if(E.is_in(polypair.second)){
            returnpair.second = polypair.second;
        }else{
            if(ellinter) //if the intersection with ball is already computed
                returnpair.second = ellpair.second;
            else returnpair.second = (E.line_intersect(r,v)).second;
        }
        return returnpair;
    }
    
    std::pair<K,K> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init
                                          ){

        std::pair<K,K> polypair = P.line_intersect_coord(r,r_prev,rand_coord,rand_coord_prev,lamdas,init);
        std::pair<K,K> ellpair = E.line_intersect_coord(r,rand_coord);
        return std::pair<K,K> (std::min(polypair.first,ellpair.first),
                                 std::max(polypair.second,ellpair.second));
    }
    
};


template <class T1 , class T2 >
class BallPolyIntersectEll {
private:
    T1 BP;
    T2 E;
    typedef typename T2::K 	K;
public:
    BallPolyIntersectEll(T1 &BPin, T2 &Ein) : BP(BPin), E(Ein) {};
    
    T1 first() { return BP; }
    T2 second() { return E; }
    
    int is_in(Point p){
        //std::cout << "calling is in"<<std::endl;
        if(BP.is_in(p)==-1)
            return E.is_in(p);
        return 0;
    }
    
    int num_of_hyperplanes(){
        return BP.num_of_hyperplanes();
    }

    int dimension(){
        return BP.dimension();
    }
    
    std::pair<Point,Point> line_intersect(Point r,
                                          Point v){

        std::pair<Point,Point> Bpolypair = BP.line_intersect(r,v);
        std::pair<Point,Point> returnpair;
        std::pair<Point,Point> ellpair;
        bool ellinter=false;

        //check the first intersection point if it is inside ball
        if(E.is_in(Bpolypair.first)){
            //std::cout<<"inside ball 1, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.first<<std::endl;
            returnpair.first = Bpolypair.first;
        }else{
            //std::cout<<"outside ball 1, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.first<<std::endl;
            ellinter=true;
            //compute the intersection with ball
            ellpair = E.line_intersect(r,v);
            returnpair.first = ellpair.first;
            //std::cout<<returnpair.first<<std::endl;
        }
        //check the second intersection point
        if(E.is_in(Bpolypair.second)){
            //std::cout<<"inside ball 2, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.second<<std::endl;
            returnpair.second = Bpolypair.second;
        }else{
            //std::cout<<"outside ball 2, radius:"<<_B.radius()<<std::endl;
            //std::cout<<polypair.second<<std::endl;
            if(ellinter) //if the intersection with ball is already computed
                returnpair.second = ellpair.second;
            else returnpair.second = (E.line_intersect(r,v)).second;
            //std::cout<<returnpair.second<<std::endl;
        }
        return returnpair;
    }

    std::pair<K,K> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init
                                          ){

        std::pair<K,K> Bpolypair = BP.line_intersect_coord(r,r_prev,rand_coord,rand_coord_prev,lamdas,init);
        std::pair<K,K> ellpair = E.line_intersect_coord(r,rand_coord);
        return std::pair<K,K> (std::min(Bpolypair.first,ellpair.first),
                                 std::max(Bpolypair.second,ellpair.second));
    }
    
    
    
};

#endif
