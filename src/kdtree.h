#ifndef KDTREE_H
#define KDTREE_H

#include "nr3.h"
#include "pointbox.h"


template<Int DIM> struct Boxnode : Box<DIM> { //kdtree.h
//Node in a binary tree of boxes containing points. See text for details.

    Int mom, dau1, dau2, ptlo, pthi;

    Boxnode() {}

    Boxnode(Point<DIM> mylo, Point<DIM> myhi,
            Int mymom, Int myd1, Int myd2, Int myptlo, Int mypthi) :
    Box<DIM>(mylo, myhi), mom(mymom), dau1(myd1),
      dau2(myd2), ptlo(myptlo), pthi(mypthi) {}
};

template<Int DIM> struct KDtree
{ //kdtree.h
//Structure for implementing a KD tree.

    static const Doub BIG; //Size of the root box, value set below.
    Int nboxes, npts; //Number of boxes, number of points.
    vector< Point<DIM> > &ptss; //Reference to the vector of points in the KD tree.
    Boxnode<DIM> *boxes; //The array of Boxnodes that form the tree.
    VecInt ptindx, rptindx; //Index of points (see text), and reverse index.
    Doub *coords; //Point coordinates rearranged contiguously.
    KDtree(vector< Point<DIM> > &pts); //Constructor.
    ~KDtree() {delete [] boxes;}

    //Next, utility functions for use after the tree is constructed. See below.

    Doub disti(Int jpt, Int kpt);
    Int locate(Point<DIM> pt);
    Int locate(Int jpt);

    //Next, applications that use the KD tree. See text.

    Int nearest(Int jpt);
    Int nearest(Point<DIM> pt);
    void nnearest(Int jpt, Int *nn, Doub *dn, Int n);
    static void sift_down(Doub *heap, Int *ndx, Int nn); ///Used by nnearest.
    Int locatenear(Point<DIM> pt, Doub r, Int *list, Int nmax);
};
template<Int DIM> const Doub KDtree<DIM>::BIG(1.0e99);

Int selecti(const Int k, Int *indx, Int n, Doub *arr)
//Permutes indx[0..n-1] to make arr[indx[0..k-1]]  arr[indx[k]]  arr[indx[k+1..n-1]].
//The array arr is not modified. See comments in the routine select.
{
    Int i,ia,ir,j,l,mid;
    Doub a;
    l=0;
    ir=n-1;
    for (;;)
    {
        if (ir <= l+1)
        {
            if (ir == l+1 && arr[indx[ir]] < arr[indx[l]])
                SWAP(indx[l],indx[ir]);
            return indx[k];
        } else
        {
            mid=(l+ir) >> 1;
            SWAP(indx[mid],indx[l+1]);
            if (arr[indx[l]] > arr[indx[ir]])
                SWAP(indx[l],indx[ir]);

            if (arr[indx[l+1]] > arr[indx[ir]])
                SWAP(indx[l+1],indx[ir]);

            if (arr[indx[l]] > arr[indx[l+1]])
                SWAP(indx[l],indx[l+1]);
            i=l+1;
            j=ir;
            ia = indx[l+1];
            a=arr[ia];
            for (;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if (j < i)
                    break;
                SWAP(indx[i],indx[j]);
            }
            indx[l+1]=indx[j];
            indx[j]=ia;
            if (j >= k)
                ir=j-1;
            if (j <= k)
                l=i;
        }
    }
}



template<Int DIM> KDtree<DIM>::KDtree(vector< Point<DIM> > &pts) : //kdtree.h
    ptss(pts), npts(pts.size()), ptindx(npts), rptindx(npts) {
    //Construct a KD tree from a vector of points.
    Int ntmp,m,k,kk,j,nowtask,jbox,np,tmom,tdim,ptlo,pthi;
    Int *hp;
    Doub *cp;
    Int taskmom[50], taskdim[50];// Enough stack for 250 points!
    for (k=0; k<npts; k++)
        ptindx[k] = k; //Initialize the index of points.
    //Calculate the number of boxes and allocate memory for them.
    m = 1;
    for (ntmp = npts; ntmp; ntmp >>= 1)
    {
        m <<= 1;
    }
    nboxes = 2*npts - (m >> 1);
    if (m < nboxes)
        nboxes = m;
    nboxes--;
    boxes = new Boxnode<DIM>[nboxes];
    //Copy the point coordinates into a contiguous array.
    coords = new Doub[DIM*npts];

    for (j=0, kk=0; j<DIM; j++, kk += npts)
    {
        for (k=0; k<npts; k++)
            coords[kk+k] = pts[k].x[j];
    }
    //Initialize the root box and put it on the task list for subdivision.
    Point<DIM> lo(-BIG,-BIG,-BIG), hi(BIG,BIG,BIG); //Syntax OK for 2-D too.
    boxes[0] = Boxnode<DIM>(lo, hi, 0, 0, 0, 0, npts-1);
    jbox = 0;
    taskmom[1] = 0;// Which box.
    taskdim[1] = 0;// Which dimension.
    nowtask = 1;
    while (nowtask) { //Main loop over pending tasks.
        tmom = taskmom[nowtask];
        tdim = taskdim[nowtask--];
        ptlo = boxes[tmom].ptlo;
        pthi = boxes[tmom].pthi;
        hp = &ptindx[ptlo]; //Points to left end of subdivision.
        cp = &coords[tdim*npts]; //Points to coordinate list for current dim.
        np = pthi - ptlo + 1;// Number of points in the subdivision.
        kk = (np-1)/2; //Index of last point on left (boundary point).
        (void) selecti(kk,hp,np,cp); //Here is where all the work is done.
        //Now create the daughters and push them onto the task list if they need further subdividing.
        hi = boxes[tmom].hi;
        lo = boxes[tmom].lo;
        hi.x[tdim] = lo.x[tdim] = coords[tdim*npts + hp[kk]];
        boxes[++jbox] = Boxnode<DIM>(boxes[tmom].lo,hi,tmom,0,0,ptlo,ptlo+kk);
        boxes[++jbox] = Boxnode<DIM>(lo,boxes[tmom].hi,tmom,0,0,ptlo+kk+1,pthi);
        boxes[tmom].dau1 = jbox-1;
        boxes[tmom].dau2 = jbox;

        if (kk > 1)
        {
            taskmom[++nowtask] = jbox-1;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
        if (np - kk > 3)
        {
            taskmom[++nowtask] = jbox;
            taskdim[nowtask] = (tdim+1) % DIM;
        }
    }
    for (j=0; j<npts; j++)
        rptindx[ptindx[j]] = j; //Create reverse index.
    delete [] coords; //Don’t need them anymore.

}

template<Int DIM> Doub KDtree<DIM>::disti(Int jpt, Int kpt)
{
    //Returns the distance between two points in the kdtree given their indices in the array of points,
    //but returns a large value if the points are identical.
    if (jpt == kpt)
        return BIG;
    else
        return dist(ptss[jpt], ptss[kpt]);
}

template<Int DIM> Int KDtree<DIM>::locate(Point<DIM> pt)
{
    //Given an arbitrary point pt, return the index of which kdtree box it is in.
    Int nb,d1,jdim;
    nb = jdim = 0; //Start with the root box.
    while (boxes[nb].dau1)
    { //As far as possible down the tree.
        d1 = boxes[nb].dau1;
        if (pt.x[jdim] <= boxes[d1].hi.x[jdim])
            nb=d1;
        else
            nb=boxes[nb].dau2;
        jdim = ++jdim % DIM; //Increment the dimension cyclically.
    }
    return nb;
}

template<Int DIM> Int KDtree<DIM>::locate(Int jpt)
{ //kdtree.h
    //Given the index of a point in the kdtree, return the index of which box it is in.
    Int nb,d1,jh;
    jh = rptindx[jpt]; //The reverse index tells where the point lies in the
    nb = 0; //index of points.
    while (boxes[nb].dau1)
    {
        d1 = boxes[nb].dau1;
        if (jh <= boxes[d1].pthi)
            nb=d1;
        else
            nb = boxes[nb].dau2;
    }
    return nb;

}



template<Int DIM> Int KDtree<DIM>::nearest(Point<DIM> pt) { //kdtree.h
    //Given an arbitrary location pt, return the index of the nearest point in the kdtree.
    Int i,k,nrst,ntask;
    Int task[50]; //Stack for boxes waiting to be opened.
    Doub dnrst = BIG, d;
    //First stage, we find the nearest kdtree point in same box as pt.

    k = locate(pt);// Which box is pt in?
    for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
    { //Find nearest.
        d = dist(ptss[ptindx[i]],pt);
        if (d < dnrst)
        {
            nrst = ptindx[i];
            dnrst = d;
        }
    }

    //Second stage, we traverse the tree opening only possibly better boxes.
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (dist(boxes[k],pt) < dnrst)
        { //Distance to closest point in box.
            if (boxes[k].dau1)
            { //If not an end node, put on task list.
                task[++ntask] = boxes[k].dau1;
                task[++ntask] = boxes[k].dau2;
            }
            else
            {// Check the 1 or 2 points in the box.
                for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
                {
                    d = dist(ptss[ptindx[i]],pt);
                    if (d < dnrst)
                    {
                        nrst = ptindx[i];
                        dnrst = d;
                    }
                }
            }
        }
    }
    return nrst;
}

template<Int DIM> void KDtree<DIM>::nnearest(Int jpt, Int *nn, Doub *dn, Int n)
//Given the index jpt of a point in a kdtree, return a list nn[0..n-1] of the indices of the n
//points in the tree nearest to point j, and a list dd[0..n-1] of their distances.
{
    Int i,k,ntask,kp;
    Int task[50]; //Stack for boxes to be opened.
    Doub d;
    if (n > npts-1)
        throw("too many neighbors requested");
    for (i=0; i<n; i++)
        dn[i] = BIG;
    //Find smallest mother box with enough points to initialize the heap.
    kp = boxes[locate(jpt)].mom;
    while (boxes[kp].pthi - boxes[kp].ptlo < n)
        kp = boxes[kp].mom;
    //  Examine its points and save the n closest.
    for (i=boxes[kp].ptlo; i<=boxes[kp].pthi; i++)
    {
        if (jpt == ptindx[i])
            continue;
        d = disti(ptindx[i],jpt);
        if (d < dn[0])
        {
            dn[0] = d;
            nn[0] = ptindx[i];
            if (n>1)
                sift_down(dn,nn,n); //Maintain the heap structure.
        }
    }
    //Now we traverse the tree opening only possibly better boxes.
    task[1] = 0;
    ntask = 1;

    while (ntask)
    {
        k = task[ntask--];
        if (k == kp)
            continue; //Don’t redo the box used to initialize.
        if (dist(boxes[k],ptss[jpt]) < dn[0])
        {
            if (boxes[k].dau1)
            { //If not an end node, put on task list.
                task[++ntask] = boxes[k].dau1;
                task[++ntask] = boxes[k].dau2;
            } else
            { //Check the 1 or 2 points in the box.
                for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
                {
                    d = disti(ptindx[i],jpt);
                    if (d < dn[0])
                    {
                        dn[0] = d;
                        nn[0] = ptindx[i];
                        if (n>1)
                            sift_down(dn,nn,n); //Maintain the heap.
                    }
                }
            }
        }
    }
    return;

}

template<Int DIM> void KDtree<DIM>::sift_down(Doub *heap, Int *ndx, Int nn) { //kdtree.h
    //Fix heap[0..nn-1], whose first element (only) may be wrongly filed. Make a corresponding
    //permutation in ndx[0..nn-1]. The algorithm is identical to that used by sift_down in hpsort.
    Int n = nn - 1;
    Int j,jold,ia;
    Doub a;
    a = heap[0];
    ia = ndx[0];
    jold = 0;
    j = 1;
    while (j <= n)
    {
        if (j < n && heap[j] < heap[j+1])
            j++;
        if (a >= heap[j])
            break;
        heap[jold] = heap[j];
        ndx[jold] = ndx[j];
        jold = j;
        j = 2*j + 1;
    }
    heap[jold] = a;
    ndx[jold] = ia;

}

template<Int DIM> Int KDtree<DIM>::locatenear(Point<DIM> pt, Doub r, Int *list, Int nmax)
{
    //Given a point pt and radius r, returns a value nret such that list[0..nret-1] is a list of all
    //kdtree points within a radius r of pt, up to a user-specified maximum of nmax points.
    Int k,i,nb,nbold,nret,ntask,jdim,d1,d2;
    Int task[50];
    nb = jdim = nret = 0;
    if (r < 0.0)
        throw("radius must be nonnegative");

    //Find the smallest box that contains the ”ball” of radius r.
    while (boxes[nb].dau1)
    {
        nbold = nb;
        d1 = boxes[nb].dau1;
        d2 = boxes[nb].dau2;
        //Only need to check the dimension that divides the daughters.
        if (pt.x[jdim] + r <= boxes[d1].hi.x[jdim])
            nb = d1;
        else if (pt.x[jdim] - r >= boxes[d2].lo.x[jdim])
            nb = d2;

        jdim = ++jdim % DIM;
        if (nb == nbold)
            break;// Neither daughter encloses the ball.
    }
    //Now traverse the tree below the starting box only as needed.
    task[1] = nb;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (dist(boxes[k],pt) > r)
            continue; //Box and ball are disjoint.

        if (boxes[k].dau1) {// Expand box further when possible.
            task[++ntask] = boxes[k].dau1;
            task[++ntask] = boxes[k].dau2;
        }
        else
        { //Otherwise process points in the box.
            for (i=boxes[k].ptlo; i<=boxes[k].pthi; i++)
            {
                if (dist(ptss[ptindx[i]],pt) <= r && nret < nmax)
                    list[nret++] = ptindx[i];
                if (nret == nmax)
                    return nmax; //Not enough space!
            }
        }
    }
    return nret;
}
#endif // KDTREE_H
