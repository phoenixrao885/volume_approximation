

template <class T>
NT Vpoly_volume (T &P, vars var) {
    int n = var.n;
    NT vol;
    std::pair<Point, NT> res = compute_minball(P);

    Point xc = res.first;
    std::cout<<"center = ";
    xc.print();
    std::cout<<"radius = "<<res.second<<std::endl;

    vol = (std::pow(M_PI,n/2.0)*(res.second, n) )  / (tgamma(n/2.0+1));

    return vol;
}