/**********************************************************/

/**********************************************************/
#include "Kerf.h"

#include <math.h>
#include <stdlib.h>

#include <fstream>
#include <string>

//double kerf = 0;
//GCodeARRAY_STRUCT GfileFloatKerf[80000]; // 有割缝的浮点型切割代码
//GCodeARRAY_STRUCT GfileFloatNoKerf[80000]; // 没有割缝的浮点型切割代码
//GraphyLimit graphylimitxy;

const double INFINITESIMAL = 0.00001;//0.0000001;
const double PI = 3.14159265358979323846264338327950288419716939937510;
const double _1p4PI = 0.78539816;
const double _1p2PI = 1.57079633;
const double _3p4PI = 2.35619449;
const double _5p4PI = 3.92699081;
const double _3p2PI = 4.71238898;
const double _7p4PI = 5.49778714;
const double _2PI = 6.28318531;
const double EP = 1e-5;

unsigned char IsZero(double x) {
  return fabs(x) < INFINITESIMAL;
}

double fpmin(double x, double y) {
  return x < y ? x : y;
}


double fpmax(double x, double y) {
  return x > y ? x : y;
}

unsigned char IsEqual(double x, double y) {
  return fabs(x - y) <= fpmax(INFINITESIMAL * fpmin(fabs(x), fabs(y)), INFINITESIMAL);
}

// return x > y
unsigned char IsGreater(double x, double y) {
  return x > y + fpmax(INFINITESIMAL * fpmin(fabs(x), fabs(y)), INFINITESIMAL);
}

// return x < y
unsigned char IsLesser(double x, double y) {
  return x < y - fpmax(INFINITESIMAL * fpmin(fabs(x), fabs(y)), INFINITESIMAL);
} 

double myatan2(double x, double y) {
  double ang = 0.;
  ang = (y != 0 || x != 0) ? atan2(x, y) : 0;
  if (ang < 0) {
    ang += _2PI;
  }
  return ang;
}

double maxfabs(double a, double b) {
  return fabs(a) > fabs(b) ? fabs(a) : fabs(b);
}

double minfabs(double a, double b) {
  return fabs(a) > fabs(b) ? fabs(b) : fabs(a);
}

double dist(double x, double y, double x1, double y1) {
  return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
}

double sqr(double a) {
  return pow(a, 2);
}

int circle_intersect(Circle_t A, Circle_t B, Point_t *ia, Point_t *ib) {
  double a, b, c, d, k, aa, bb, cc, dd, drt;
  if (IsEqual(A.x, B.x) && IsEqual(A.y, B.y)) {
    return 5; // 同心圆
  }
  dd = dist(A.x, A.y, B.x, B.y);
  if (A.r + B.r + EP < dd) {
    return 1;  // 相离
  }
  a = B.x - A.x;
  b = B.y - A.y;
  c = B.r;
  k = A.r;
  d = sqr(c) - sqr(k) - sqr(a) - sqr(b);

  aa = 4 * sqr(a) + 4 * sqr(b);
  bb = 4 * b * d;
  cc = sqr(d) - 4 * sqr(a) *sqr(k);

  drt = sqr(bb) - 4 * aa * cc;
  if (drt < 0) {
    return 5; // 无解
  }
  drt = sqrt(drt);
  ia->y = (-bb + drt) / 2 / aa;
  ib->y = (-bb - drt) / 2 / aa;
  if (fabs(a) < EP) {
    ia->x = sqrt(sqr(k) - sqr(ia->y));
    ib->x = -ia->x;
  } else {
    ia->x = (2 * b * ia->y + d) / -2 / a;
    ib->x = (2 * b * ib->y + d) / -2 / a;
  }
  ia->x += A.x;
  ia->y += A.y;
  ib->x += A.x;
  ib->y += A.y;
  if (fabs(ia->y - ib->y) < EP) {
    if (fabs(A.r + B.r - dd) < EP) {
      return 2; // 外切一个交点
    }
    if (fabs(dd - (maxfabs(A.r, B.r) - minfabs(A.r, B.r))) < EP) {
      return 3; // 内切一个交点
    }
  }
  return 4; // 2个交点
}

/**********************************************************/

/**********************************************************/
int CalcRoots(double a, double b, double c, double *root1, double *root2) {
  if (IsZero(a)) {
    if (IsZero(b)) {
      return IsZero(c) ? -1 : 0;
    } else {
      *root1 = -c / b;
      return 1;
    }
  } else {
    double delta = b * b - 4 * a * c;
    
    if (IsLesser(delta, 0.)) { // no root
      return 0;
    } else if (IsZero(delta)) {
      *root1 = -b / (2 * a);
      return 1;
    } else {
      delta = sqrt(delta);
      *root1 = (-b + delta) / (2 * a);
      *root2 = (-b - delta) / (2 * a);
      return 2;
    }
  }
}
//lead_type 0--交点靠近圆弧引入点 1--交点靠近圆弧引出点
/**********************************************************/
int if_in_arc_ccw(double ang_temp,double ang_str,double ang_end,int lead_type)
{
	/*ang_temp用来看是否在圆弧内，已经将角度变为+但360度附近还是*/
	/*要注意，参考画圆部分*/
	/*如果str>end则temp>str||temp<end即可*/
	/*如果str<end则str<temp<end即可*/
	/*
  if (ang_str<0) ang_str+=6.28;
	if (ang_end<0) ang_end+=6.28;
	if (ang_temp<0) ang_temp+=6.28;
	if (ang_str>ang_end)
	{
		if (ang_temp>ang_str||ang_temp<ang_end)  return 0;
	}
	else if (ang_str<ang_end)
	{
		if (ang_temp>=ang_str&&ang_temp<=ang_end)  return 0;
	}
	else
		return 0;
	
	return 1;
  */
	if (IsLesser(ang_str ,0.)) {
    ang_str += 2. * PI;
  }
  if (IsLesser(ang_end ,0.)) {
    ang_end += 2. * PI;
  }
  if (IsLesser(ang_temp, 0.)) {
    ang_temp += 2. * PI;
  }
  if (IsEqual(ang_str, ang_end)) {
    ang_end += 2. * PI;
    if(ang_temp>ang_str && lead_type==1)
      ang_temp += 2. * PI;
  }
  if (IsGreater(ang_str, ang_end)) {
    if (IsGreater(ang_temp, ang_str) || IsLesser(ang_temp, ang_end)) {
      return 0;
    }
  } else if (IsLesser(ang_str, ang_end)) {
    if (!IsLesser(ang_temp, ang_str) && !IsGreater(ang_temp, ang_end)) {
      return 0;
    }
  } else {
    return 0;
  }

  return 1;
}
//lead_type 0--交点靠近圆弧引入点 1--交点靠近圆弧引出点
int if_in_arc_cw(double ang_temp,double ang_str,double ang_end,int lead_type)
{
	/*ang_temp用来看是否在圆弧内，已经将角度变为+但360度附近还是*/
	/*要注意，参考画圆部分*/
	/*如果str>end则str>temp>end即可*/
	/*如果str<end则temp<str||temp>end即可*/
	/*
	if (ang_str<0) ang_str+=_2PI;
	if (ang_end<0) ang_end+=_2PI;
	if (ang_temp<0) ang_temp+=_2PI;
	
	if (ang_str>ang_end)
	{
		if (ang_temp<=ang_str&&ang_temp>=ang_end)  return 0;
	}
	else if (ang_str<ang_end)
	{
		if (ang_temp<ang_str||ang_temp>ang_end)  return 0;
	}
	else
		return 0;
	
	return 1;*/
	if (IsLesser(ang_str, 0.)) {
    ang_str += 2. * PI;
  }
  if (IsLesser(ang_end, 0.)) {
    ang_end += 2. * PI;
  }
  if (IsLesser(ang_temp, 0.)) {
    ang_temp += 2. * PI;
  }
  if (IsEqual(ang_str, ang_end)) {
    ang_end -= 2. * PI;
    if(ang_temp<ang_str&& lead_type==1)
      ang_temp -= 2. * PI;
  }
  if (IsGreater(ang_str, ang_end)) {
    if (!IsGreater(ang_temp, ang_str) && !IsLesser(ang_temp, ang_end)) {
      return 0;
    }
  } else if (IsLesser(ang_str, ang_end)) {
    if (IsLesser(ang_temp, ang_str) || IsGreater(ang_temp, ang_end)) {
      return 0;
    }
  } else {
    return 0;
  }

  return 1;
}
/*
	x0,y0 start(x,y)
	x1,y1 end(x,y)
	cx,cy end(x,y)
	px,py test(x,y)
*/
int if_point_in_arc(double x0,double y0, double x1,double y1, double cx,double cy,double px,double py)
{
	double ang_str,ang_end,ang_pt;
	double d_t1,d_t2,d_t3,d_t4;
	d_t1 = x1-cx;     // 由终点相对起点坐标转换成终点相对圆心坐标
  d_t2 = y1-cy;
  d_t3 = x0-cx;   // 由圆心相对起点坐标转换成起点相对圆心坐标
  d_t4 = y0-cy;
	ang_end = myatan2(d_t2,d_t1);  
	ang_str = myatan2(d_t4, d_t3);
	d_t3 = px-cx;   
  d_t4 = py-cy;
	ang_pt = myatan2(d_t4, d_t3);

    while(ang_str<0)
    ang_str+=_2PI;
    while(ang_str>_2PI)
    ang_str-=_2PI;

    while(ang_end<0)
    ang_end+=_2PI;
    while(ang_end>_2PI)
    ang_end-=_2PI;
    
    while(ang_pt<0)
    ang_pt+=_2PI;
    while(ang_pt>_2PI)
    ang_pt-=_2PI;
 	/*2014.6.16if(ang_end<ang_str) {
     ang_end += _2PI;
    }
    d_t1 = (ang_end-ang_pt)*(ang_pt-ang_str);
	if(d_t1>=0)
		return 1;
	else
		return 0;
		*/
    d_t1 = ang_end-ang_str;
    if(d_t1 <= EP)
    	d_t1 += _2PI;

   	d_t2 = ang_pt-ang_str;
    if(d_t2<=0)
    	d_t2+=_2PI;
    if(d_t1>=d_t2||fabs(ang_pt-ang_end)<=0.00002||fabs(ang_pt-ang_str)<=0.00002)
    return 1;
    else
    return 0;
	
}


/*
int if_in_line(double ang1x,double ang1y,double ang2x,double ang2y,double ang3x,double ang3y)
{
	double R,R2,R3;

	R=pow(ang2x-ang3x,2)+pow(ang2y-ang3y,2);
	R2=pow(ang1x-ang2x,2)+pow(ang1y-ang2y,2);
	R3=pow(ang1x-ang3x,2)+pow(ang1y-ang3y,2);

	if (R>R2&&R>R3)
		return 0;
	else
		return 1;
}*/

int point_in_line(double pointx,double pointy,double x0,double y0,double x,double y)
{
	double R,R2;
	R = (x-pointx)*(pointx - x0);
	R2 =(y-pointy)*(pointy-y0);
	if(R>=0&&R2>=0)  //点在线段上
		return 1;
	else
		return 0;
	
}

int if_in_line_line(double ang1x,double ang1y,double ang2x,double ang2y,double ang3x,double ang3y)
{
	double R,R2,R3;

	R=pow(ang2x-ang3x,2)+pow(ang2y-ang3y,2);
	R2=pow(ang1x-ang2x,2)+pow(ang1y-ang2y,2);
	R3=pow(ang1x-ang3x,2)+pow(ang1y-ang3y,2);

	if (R>R2&&R>R3)
		return 0;
	else if (R2>R3)
	       return 0;
	else
		return 1;
}
/*
返回值 :
	1,5-相离
	2-外切
	3-内切
	4-两个交点
*/
int CircleIntersector(double x1, double y1, 
                      double r1, double x2, 
                      double y2, double r2, 
                      double *ix1, double *iy1, 
                      double *ix2, double *iy2) 
{
  double dd = dist (x1, y1, x2, y2);
  // double r;
  // double c1;
  // double c2;
  // double s1;
  // double s2;
  double a;
  double b;
  double c;
  double p;
  double q;
  int rtn;
  
  Circle_t A, B;
  Point_t PointA, PointB;
  if (IsLesser(r1 + r2, dd)) {
    return 1; // 相离
  }

  if (IsLesser(dd, fabs(r1 - r2)) || (IsEqual(x1, x2) && IsEqual(y1, y2) && IsEqual(r1, r2))) {
    return 5; // 相含
  }

  a = 2. * r1 * (x1 - x2);
  b = 2. * r1 * (y1 - y2);
  c = sqr(r2) - sqr(r1) - sqr(x1 - x2) - sqr(y1 - y2);
  p = sqr(a) + sqr(b);
  q = -2. * a * c;
  
  if((1 - sqr(-q / p / 2))<0) // 几乎相切或相交，但实际相离的情况
		return 1;
  // 相切
  if (IsEqual(dd, fabs(r1 - r2)) || IsEqual(dd, r1 + r2)) {
    if(IsEqual(dd, fabs(r1 - r2)))
    {
        if(r1>r2)
        {
            *ix1 = *ix2 = x2+r2*(x2-x1)/(r1-r2);
            *iy1 = *iy2 = y2+r2*(y2-y1)/(r1-r2);
        }
        else
        {
            *ix1 = *ix2 = x1+r1*(x1-x2)/(r2-r1);
            *iy1 = *iy2 = y1+r1*(y1-y2)/(r2-r1);
        }
        return 3;
    }
    else
    {
        *ix1 = *ix2 = x1+r1*(x2-x1)/(r1+r2);
        *iy1 = *iy2 = y1+r1*(y2-y1)/(r1+r2);
        return 2;
    }
    
   /* c1 = -q / p / 2.;
    s1 = sqrt(1 - sqr(c1));
    *ix1 =  r1 * c1 + x1;
    *iy1 =  r1 * s1 + y1;
    if (!IsEqual(sqr(*ix1 - x2) + sqr(*iy1 - y2), sqr(r2))) {
      *iy1 = -r1 * s1 + y1;
    }
    *ix2  = *ix1;
    *iy2 = *iy1;
    return IsEqual(dd, fabs(r1 - r2)) ? 3 : 2;*/
  }
  A.x = x1;
  A.y = y1;
  A.r = r1;
  B.x = x2;
  B.y = y2;
  B.r = r2;
  rtn = circle_intersect(A,B,&PointA,&PointB);
  if(rtn ==1)
  {
    //同心圆
    return 5;
  }
  else if(rtn==4)//两个交点
  {
      *ix1 = PointA.x;
      *iy1 = PointA.y;
      *ix2 = PointB.x;
      *iy2 = PointB.y;
      return 4;
  }
  else if(rtn==2)//外切一个交点
  {
      *ix1 = PointA.x;
      *iy1 = PointA.y;
      *ix2 = PointB.x;
      *iy2 = PointB.y;
      return 2;
  }
  else if(rtn ==3 )//内切一个交点
  {
      *ix1 = PointA.x;
      *iy1 = PointA.y;
      *ix2 = PointB.x;
      *iy2 = PointB.y;
      return 3;
  }
  else
  {
      //fault
      return 1;
  
  }

}
double ang_fabs(double ang)
{
	if(ang>0)
	{
		if(ang>1.57)
			ang-=3.14;
	}
	else if (ang<0)
	{
		if (ang<-1.57)
			ang+=3.14;
	}
	return ang;
}
double mysqrt(double x)
{
	if (x<0.001) x=0;
	return sqrt(x);
}

/*
返回值 :
  0 平行
  1 有交点
*/
int TwoLineIsIntersect(double   x0,   double   y0,   double   x1,   double   y1,   double   x2,   double   y2,   double   x3,   double   y3,   double   *InterX,   double   *InterY)   
{ //两条线段是否相交X0X1   AND   X1X2   
  //double   x,   y;   
  #define Min(v0,v1) ((v0>v1) ? v1 : v0)
  #define Max(v0,v1) ((v0>v1) ? v0 : v1)
  double   Minx01   =   Min(x0,   x1);   
  double   Miny01   =   Min(y0,   y1);   
  double   Minx23   =   Min(x2,   x3);   
  double   Miny23   =   Min(y2,   y3);   
  double   Maxx01   =   Max(x0,   x1);   
  double   Maxy01   =   Max(y0,   y1);   
  double   Maxx23   =   Max(x2,   x3);   
  double   Maxy23   =   Max(y2,   y3);   
  double   k1;   
  double   k2;   
  //double   Den;
	if(fabs(x0-x1)>=fabs(y0-y1)&&fabs(x2-x3)>=fabs(y2-y3))
	{
		k1 = (y0-y1)/(x0-x1);
		k2 = (y2-y3)/(x2-x3);
		if(IsEqual(k1,k2))
		{
			if(IsEqual(x1,x2)&&IsEqual(y1,y2))
			{
				*InterX = x1;
				*InterY = y1;
				return 1;
			}
			else
				return 0;
		}
		
	}
	else if(fabs(x0-x1)<=fabs(y0-y1)&&fabs(x2-x3)<=fabs(y2-y3))
	{
		k1 = (x0-x1)/(y0-y1);
		k2 = (x2-x3)/(y2-y3);
		if(IsEqual(k1,k2))
		{
			if(IsEqual(x1,x2)&&IsEqual(y1,y2))
			{
				*InterX = x1;
				*InterY = y1;
				return 1;
			}
			else
				return 0;
		}
		
	}

	//else
	{
		if(IsZero(x1-x0))
		{
			*InterX = x1;
			if(IsZero(x3-x2))
				return 0;
			
			k1 = (y3-y2)/(x3-x2);
			*InterY = y2+k1*(*InterX-x2);
			return 1;
		}
		else if(IsZero(x3-x2))
		{
			*InterX = x2;
			if(IsZero(x1-x0))
				return 0;
			
			k1 = (y1-y0)/(x1-x0);
			*InterY = y0+k1*(*InterX-x0);
			return 1;
		}
		else
		{
			k1 = (y0-y1)/(x0-x1);
			k2 = (y2-y3)/(x2-x3);
			*InterX = (y2-y0 +k1*x0-k2*x2)/(k1-k2);
			//*InterY = y0+k1*(*InterX-x0);
			*InterY  = (k1*k2*(x2-x0)+k2*y0-k1*y2)/(k2-k1);
			return 1;
		}
	}
	return -1;
  }
/*int   TwoLineIsIntersect(double   x0,   double   y0,   double   x1,   double   y1,   double   x2,   double   y2,   double   x3,   double   y3,   double   *InterX,   double   *InterY)   
  {   //两条线段是否相交X0X1   AND   X1X2   
        double   x,   y;   
        #define Min(v0,v1) ((v0>v1) ? v1 : v0)
		#define Max(v0,v1) ((v0>v1) ? v0 : v1)
        double   Minx01   =   Min(x0,   x1);   
        double   Miny01   =   Min(y0,   y1);   
        double   Minx23   =   Min(x2,   x3);   
        double   Miny23   =   Min(y2,   y3);   
        double   Maxx01   =   Max(x0,   x1);   
        double   Maxy01   =   Max(y0,   y1);   
        double   Maxx23   =   Max(x2,   x3);   
        double   Maxy23   =   Max(y2,   y3);   
	 float   k1   =   (y1-y0)/(x1-x0);   
        float   k2   =   (y3-y2)/(x3-x2);   
        float   Den   =   (y1-y0)*(x3-x2)   -   (y3-y2)*(x1-x0);   
	if(fabs(x0-x1)>fabs(y0-y1)&&fabs(x2-x3)>fabs(y2-y3))
	{
		k1 = (y0-y1)/(x0-x1);
		k2 = (y2-y3)/(x2-x3);
		if(IsEqual(k1,k2))
			return 0;
	}
	if(fabs(x0-x1)<fabs(y0-y1)&&fabs(x2-x3)<fabs(y2-y3))
	{
		k1 = (x0-x1)/(y0-y1);
		k2 = (x2-x3)/(y2-y3);
		if(IsEqual(k1,k2))
			return 0;
	}
          
        if(!IsEqual(x1,x0)   &&  !IsEqual(x2,x3))   
        {   

                if(IsEqual(k1,k2))   
                {   //平行不相交   
                      float   d1   =   abs(y0*(x1-x0)-x0*(y1-y0)-y2*(x3-x2)+x2*(y3-y2));   //距离公式d   =   abs(c1-c2)   /   sqrt(a*a+b*b)   
                      if(d1==0)   
                      {//直线重合   
                            if((x2>Minx01   &&   x2<Maxy01   &&   y2>Miny01   &&   y2<Maxy01)   ||   (x3>Minx01   &&   x3<Maxy01   &&   y3>Miny01   &&   y3<Maxy01)   
                            ||   (x0>Minx23   &&   x0<Maxy23   &&   y0>Miny23   &&   y0<Maxy23)   ||   (x1>Minx23   &&   x1<Maxy23   &&   y1>Miny23   &&   y1<Maxy23))   
                            {     //实际碰撞问题线段重合认为相交了   
                                  return   0;//;true;   
                            }   
                            else   
                            {   
                                  return   0;//false;   
                            }   
                      }   
                      else   
                      {   
                            return   0;//false;   
                      }         
                }   
                x   =   ((y2-y0)*(x1-x0)*(x3-x2)+(y1-y0)*(x3-x2)*x0-(y3-y2)*(x1-x0)*x2)/Den;   
                y   =   ((y1-y0)*(x-x0))/(x1-x0)   +   y0;   
    
                if(Minx01<=x   &&   x<=Maxx01   &&   Miny01<=y   &&   y<=Maxy01   &&   Minx23<=x   &&   x<=Maxx23   &&   Miny23<=y   &&   y<=Maxy23)   
                {   
                      *InterX   =   x;   
                      *InterY   =   y;   
                      return   1;//true;   
                }   
        }   
        else   if(IsEqual(x1,x0)   &&  !IsEqual(x2,x3))   
        {   
                x   =   x0;   
                y   =   ((y3-y2)*(x0-x2))/(x3-x2)   +   y2;   
                if(Minx01<=x   &&   x<=Maxx01   &&   Miny01<=y   &&   y<=Maxy01   &&   Minx23<=x   &&   x<=Maxx23   &&   Miny23<=y   &&   y<=Maxy23)   
                {   
                      *InterX   =   x;   
                      *InterY   =   y;   
                      return   1;//true;   
                }   
        }   
        else   if(!IsEqual(x1,x0)   &&   IsEqual(x2,x3))   
        {   
                x   =   x2;   
                y   =   ((y1-y0)*(x2-x0))/(x1-x0)   +   y0;   
                if(Minx01<=x   &&   x<=Maxx01   &&   Miny01<=y   &&   y<=Maxy01   &&   Minx23<=x   &&   x<=Maxx23   &&   Miny23<=y   &&   y<=Maxy23)   
                {   
                      *InterX   =   x;   
                      *InterY   =   y;   
                      return  1;//;   
                }                 
        }   
        return   -1;//false;   
  }*/
/*
	返回值:
		0-没有交点，即相离
		1-一个交点，即相切
		2-二两个交点，即相交
		
*/
int LineCircleIntersect(double x0,double y0,double x,double y, double xc0,double yc0,double Radius,double *InterX1,double *InterY1,double *InterX2,double *InterY2)
{
	int res=0;
	double k1=1,temp1;
	double x_intersect1 = 0,y_intersect1 = 0;
	double x_intersect2 = 0,y_intersect2 = 0;
	double a,b,c,bb4ac;//ax*x+bx+c=0
	/*x=(-b+-mysqrt(b*b-4ac))/2a*/
	/*先判断b*b-4ac>0?*/
	/*>0:直线与园相交*/
	/*=0:直线与园相切*/
	/*<0:直线与园相离即外割逢*/

	if(fabs(x-x0)>fabs(y-y0))
	{
		k1 = (y-y0)/(x-x0);
		//y = k1*x-k1*x0+y0代入圆方程,求解X
		temp1 = y0 - k1*x0;
		a = 1+k1*k1;
		b = -2*xc0+2*k1*(temp1-yc0);
		c = pow(xc0,2)  + pow(temp1-yc0,2)-pow(Radius,2);
              bb4ac = pow(b,2)-4*a*c;
	       if(bb4ac>0.1) //相交
		{
			res = 2;
			*InterX1 = (-b+sqrt(bb4ac))/2/a;
			*InterX2 = (-b-sqrt(bb4ac))/2/a;
			*InterY1 = k1*(*InterX1)+temp1;
			*InterY2 = k1*(*InterX2)+temp1;
		}
		else if(bb4ac<0.001)  //相离
		{
			res = 0;
		}
		else 
		{
			res = 1;
			*InterX1 = *InterX2 = (-b)/2/a;
			*InterY1 = *InterY2 = k1*(*InterX2)+temp1;
		}
	}
	else
	{
		k1 = (x-x0)/(y-y0);
		//y = k1*x-k1*x0+y0代入圆方程,求解X
		temp1 = x0 - k1*y0;
		a = 1+k1*k1;
		b = -2*yc0+2*k1*(temp1-xc0);
		c = pow(yc0,2)  + pow(temp1-xc0,2)-pow(Radius,2);
              bb4ac = pow(b,2)-4*a*c;
	       if(bb4ac>0.1) //相交
		{
			res = 2;
			*InterY1= (-b+sqrt(bb4ac))/2/a;
			*InterY2 = (-b-sqrt(bb4ac))/2/a;
			*InterX1 = k1*(*InterY1)+temp1;
			*InterX2 = k1*(*InterY2)+temp1;
		}
		else if(bb4ac<0.001)  //相离
		{
			res = 0;
		}
		else 
		{
			res = 1;
			*InterY1 = *InterY2 = (-b)/2/a;
			*InterX1 = *InterX2 = k1*(*InterY2)+temp1;
		}
	}
	return res;

}

Kerf::Kerf(): kerf(0.), is_absolute(false) {}

Kerf::~Kerf() {
}

void Kerf::SetKerfValue(double kerf_value) {
  kerf = kerf_value;
}

//GCodeARRAY_STRUCT GfileFloatTemp[GfileTotalRows];
void Kerf::CircleCheFen(std::vector<GCodeARRAY_STRUCT> &GCodeArry)
{
	double StartAngle,EndAngle,StepAngle,CurrentAngle;
	double d_t1,d_t2,d_t3,d_t4;
  std::vector<GCodeARRAY_STRUCT>::iterator GCodeArryPtrSrc = GCodeArry.begin();
	std::vector<GCodeARRAY_STRUCT>::iterator GCodeArryCheFen = GfileFloatKerf.begin();
  //GfileFloatTemp;
	
	while(GCodeArryPtrSrc->Name!=M02)
	{
		*GCodeArryCheFen = *GCodeArryPtrSrc;
		if(GCodeArryPtrSrc->Name==G02||GCodeArryPtrSrc->Name==G03)
		{
			if(GCodeArryPtrSrc->R>50000)
			{
				d_t1 = GCodeArryPtrSrc->X-GCodeArryPtrSrc->I;     // 由终点相对起点坐标转换成终点相对圆心坐标
  			d_t2 = GCodeArryPtrSrc->Y-GCodeArryPtrSrc->J;
  			d_t3 = GCodeArryPtrSrc->X0 -GCodeArryPtrSrc->I;   // 由圆心相对起点坐标转换成起点相对圆心坐标
    		d_t4 = GCodeArryPtrSrc->Y0 -GCodeArryPtrSrc->J;
    		StartAngle= myatan2(d_t4, d_t3);
    		EndAngle= myatan2(d_t2, d_t1);
	    		
  		  //StepAngle = 0.00087-0.00076*((GCodeArryPtrSrc->R-50000)/100000000);
  		  if(GCodeArryPtrSrc->R>500000)
  		  {
  			  //StepAngle = 5/GCodeArryPtrSrc->R;// 1-0.1 2-0.2
  			  StepAngle = 2*sqrt(2*0.1*GCodeArryPtrSrc->R-0.1*0.1)/GCodeArryPtrSrc->R;// pricise 0.1mm
  		  }
  		  else
  			  StepAngle = 2*((2*sqrt((100*GCodeArryPtrSrc->R)*(GCodeArryPtrSrc->R)-(GCodeArryPtrSrc->R*10-1)*(GCodeArryPtrSrc->R*10-1)))/GCodeArryPtrSrc->R/10);// 1-0.1 2-0.2
  		  //StepAngle = arccos(1-0.1/GCodeArryPtrSrc->R);
  		  if(GCodeArryPtrSrc->Name==G02)
  		  {
  			  if(StartAngle<=EndAngle+EP)
  				  StartAngle+=_2PI;
  		  }
  		  else if(GCodeArryPtrSrc->Name==G03)
  		  {
  			  if(StartAngle>=EndAngle+EP)
  				  StartAngle-=_2PI;
  		  }
  		  CurrentAngle = StartAngle;
  		  do{
  			  if(GCodeArryPtrSrc==GCodeArry.begin())
  			  {
  				  GCodeArryCheFen->X0=0;
  				  GCodeArryCheFen->Y0=0;
  			  }
  			  else
  			  *GCodeArryCheFen=*(GCodeArryCheFen-1);
  			  GCodeArryCheFen->Name = G01;
  			  GCodeArryCheFen->ShowLine = GCodeArryPtrSrc->ShowLine;
  			  //if(GCodeArryPtrSrc->F < 0.83)
  			  GCodeArryCheFen->F=0;
  			  GCodeArryCheFen->X0 = (GCodeArryCheFen-1)->X;
  			  GCodeArryCheFen->Y0 = (GCodeArryCheFen-1)->Y;
  			  if(GCodeArryPtrSrc->Name==G02)
  			  {
  				  d_t2 = CurrentAngle-EndAngle;
  				  if(d_t2>=_2PI)
  					  d_t2-=_2PI;
  				  else if(d_t2<=-_2PI)
  					  d_t2+=_2PI;
  				  if(d_t2<StepAngle)
  				  {
  					  GCodeArryCheFen->X = GCodeArryPtrSrc->X;
  					  GCodeArryCheFen->Y = GCodeArryPtrSrc->Y;
  					  CurrentAngle = EndAngle;
  				  }
  				  else
  				  {
  					  CurrentAngle -= StepAngle;
    				  if(CurrentAngle<=-_2PI)
    					  CurrentAngle-=_2PI;
    				  d_t2 = fabs(CurrentAngle-EndAngle);
    				  if(d_t2>=_2PI)
    					  d_t2-=_2PI;
    				  else if(d_t2<=-_2PI)
    					  d_t2+=_2PI;
    				  if(d_t2<=StepAngle)
    				  {
    					  GCodeArryCheFen->X = GCodeArryPtrSrc->X;
    					  GCodeArryCheFen->Y = GCodeArryPtrSrc->Y;
    					  CurrentAngle = EndAngle;
    				  }
    				  else
    				  {
    					  GCodeArryCheFen->X = GCodeArryPtrSrc->I+GCodeArryPtrSrc->R*cos(CurrentAngle);
    					  GCodeArryCheFen->Y = GCodeArryPtrSrc->J+GCodeArryPtrSrc->R*sin(CurrentAngle);
    				  }
  				  }
    				
  			  }
  			  else
  			  {
  				  d_t2 = EndAngle-CurrentAngle;
  				  if(d_t2>=_2PI)
  					  d_t2-=_2PI;
  				  else if(d_t2<=-_2PI)
  					  d_t2+=_2PI;
  				  if(d_t2<StepAngle)
  				  {
  					  GCodeArryCheFen->X = GCodeArryPtrSrc->X;
  					  GCodeArryCheFen->Y = GCodeArryPtrSrc->Y;
  					  CurrentAngle = EndAngle;
  				  }
  				  else
  				  {
    				  CurrentAngle += StepAngle;
    				  if(CurrentAngle>=_2PI)
    				  CurrentAngle-=_2PI;
    				  d_t2 = fabs(CurrentAngle-EndAngle);
    				  if(d_t2>=_2PI)
    					  d_t2-=_2PI;
    				  else if(d_t2<=-_2PI)
    					  d_t2+=_2PI;
    				  if(d_t2<=StepAngle)
    				  {
    					  GCodeArryCheFen->X = GCodeArryPtrSrc->X;
    					  GCodeArryCheFen->Y = GCodeArryPtrSrc->Y;
    					  CurrentAngle = EndAngle;
    				  }
    				  else
    				  {
    					  GCodeArryCheFen->X = GCodeArryPtrSrc->I + GCodeArryPtrSrc->R*cos(CurrentAngle);
    					  GCodeArryCheFen->Y = GCodeArryPtrSrc->J + GCodeArryPtrSrc->R*sin(CurrentAngle);
    				  }
  				  }
  			  }
  			  GCodeArryCheFen++;
  		  } while(CurrentAngle != EndAngle); 
		    GCodeArryCheFen--;
			}
		}
		GCodeArryPtrSrc++;
		GCodeArryCheFen++;
		
	}
	GCodeArryCheFen->Name = M02;
	GCodeArryPtrSrc = GCodeArry.begin();
	GCodeArryCheFen = GfileFloatKerf.begin();//GfileFloatTemp;
	while(GCodeArryCheFen->Name!=M02)
	{
		*GCodeArryPtrSrc++=*GCodeArryCheFen++;
		
	}
	GCodeArryPtrSrc->Name = M02;
	
}

int Kerf::IgnoreLittleLine(std::vector<GCodeARRAY_STRUCT> &GCodeArry) {
	std::vector<GCodeARRAY_STRUCT>::iterator GCodeArryPtrSrc = GCodeArry.begin();
	std::vector<GCodeARRAY_STRUCT>::iterator GCodeArryLittelLine = GfileFloatKerf.begin();
  //GfileFloatTemp;

	int res=0; 
  double AccLength=0;

	while(GCodeArryPtrSrc->Name!=M02)
	{
		if(GCodeArryPtrSrc->Length<0.1&&AccLength<0.5&&(GCodeArryPtrSrc->X0!=0||GCodeArryPtrSrc->Y0!=0)&&(GCodeArryPtrSrc->Name>=G00&&GCodeArryPtrSrc->Name<=G03)&&GCodeArryPtrSrc>GCodeArry.begin()&&GCodeArryPtrSrc->Name==G01&&(GCodeArryPtrSrc-1)->Name==G01)
		{
		    AccLength += GCodeArryPtrSrc->Length;

			GCodeArryLittelLine--;
			
			GCodeArryLittelLine->X=(GCodeArryPtrSrc)->X;
			GCodeArryLittelLine->Y=(GCodeArryPtrSrc)->Y;
			//GCodeArryPtrSrc++;
			res = 1;
			
		}
		else
		{
		  AccLength = 0;
			*GCodeArryLittelLine=*GCodeArryPtrSrc;
			
		}
		GCodeArryLittelLine++;
		GCodeArryPtrSrc++;
	}
	GCodeArryLittelLine->Name=M02;
	
	GCodeArryLittelLine = GfileFloatKerf.begin();//GfileFloatTemp;
	GCodeArryPtrSrc = GCodeArry.begin();
	while(GCodeArryLittelLine->Name!=M02)
	{
		*GCodeArryPtrSrc++=*GCodeArryLittelLine++;
	}
	GCodeArryPtrSrc->Name=M02;
	return res;
}

/*计算两个直线的交点坐标
返回值:
	0-没有交点
	1- 有交点

int GetLineLineIntersect(Line_t *pPreviousLine,Line_t *pNextLine,point_t *pJointPoint)
{
	int res = 0;  //无交点
	double k1,k2,Tan1,Tan2;
	double x_intersect,y_intersect;
	Tan1=myatan2(pPreviousLine.Y-pPreviousLine.X,pPreviousLine.Y0-pPreviousLine.X0);//第一条线的角度
	Tan2=myatan2(pNextLine.Y-pNextLine.X,pNextLine.Y0-pNextLine.X0);//第二条线的角度
	if (fabs(tan(Tan1)-tan(Tan2))<0.001)  break;  //平行线没有交点
	else
	{	
		k1 = (pPreviousLine->Y-pPreviousLine->Y0)/(pPreviousLine->X-pPreviousLine->X0);
		k2 = (pNextLine->Y-pNextLine->Y0)/(pNextLine->X-pNextLine->X0);
		pJointPoint->x = (k1*pPreviousLine->X0-k2*pNextLine->X0+pNextLine->Y0-pPreviousLine->Y0)/(k1-k2);
		pJointPoint->y = k1*(pJointPoint->x-pPreviousLine->X0)+pPreviousLine->Y0;
		res = 1;  //

	}
	return res;
}*/
double Kerf::GetRadius(GCodeARRAY_STRUCT &pG)
{
	return sqrt(pow(pG.X-pG.I,2)+pow(pG.Y-pG.J,2));
}
/*
	StartOrEnd:
	 0- startpoint
	 1 -endpoint
*/
double Kerf::GetTangent(GCodeARRAY_STRUCT *pG, int StartOrEnd)
{
	double Tangent;
	if(pG->Name==G01)
	{
		Tangent = myatan2(pG->Y-pG->Y0, pG->X-pG->X0);
	}
	else if(pG->Name==G02)
	{
		if(StartOrEnd==0)
			Tangent = myatan2(pG->Y0-pG->J, pG->X0-pG->I)-_1p2PI;
		else
			Tangent = myatan2(pG->Y-pG->J, pG->X-pG->I)-_1p2PI;
	}
	else if(pG->Name==G03)
	{
		if(StartOrEnd==0)
			Tangent = myatan2(pG->Y0-pG->J, pG->X0-pG->I)+_1p2PI;
		else
			Tangent = myatan2(pG->Y-pG->J, pG->X-pG->I)+_1p2PI;
	}
	return Tangent;
}
/*
    对未有增加割缝补偿的G 代码未补偿后的G代码
    只对G01, G02,G03进行计算
    kerfvalue < 0为左补偿，
    kerfvalue>0为右补偿
    返回值 的Name==M02为错误补偿 
*/
int Kerf::GetAddKerfGCode(GCodeARRAY_STRUCT &pNoKerfG, GCodeARRAY_STRUCT *AddKerfGCode,double kerfvalue,int dir)
{
	//GCodeARRAY_STRUCT KerfGCode;
	double KerfX=0,KerfY=0;
	double LastTangent, NextTangent;
	AddKerfGCode->Name = M02;
	if(fabs(kerfvalue)<EP)
		return 1;
	
	if(pNoKerfG.Name==G01)
	{
		LastTangent = myatan2(pNoKerfG.Y-pNoKerfG.Y0, pNoKerfG.X-pNoKerfG.X0);
		if(dir==G41)
		{
			KerfX = (kerfvalue)*cos(LastTangent+_1p2PI);
			KerfY = (kerfvalue)*sin(LastTangent+_1p2PI);
		}
		else
		{
			KerfX = (kerfvalue)*cos(LastTangent -_1p2PI);
			KerfY = (kerfvalue)*sin(LastTangent - _1p2PI);
		}
		*AddKerfGCode = pNoKerfG;
		AddKerfGCode->X0 += KerfX;
		AddKerfGCode->Y0 += KerfY;
		AddKerfGCode->X  += KerfX;
		AddKerfGCode->Y  += KerfY;
		
	}
	else if(pNoKerfG.Name==G02)
	{
		LastTangent = myatan2(pNoKerfG.Y0-pNoKerfG.J, pNoKerfG.X0-pNoKerfG.I)-_1p2PI;
		NextTangent = myatan2(pNoKerfG.Y-pNoKerfG.J, pNoKerfG.X-pNoKerfG.I)-_1p2PI;
		
		if(dir==G41)
		{
			KerfX = (kerfvalue)*cos(LastTangent+_1p2PI);
			KerfY = (kerfvalue)*sin(LastTangent+_1p2PI);	
		}
		else
		{
			KerfX = (kerfvalue)*cos(LastTangent-_1p2PI);
			KerfY = (kerfvalue)*sin(LastTangent-_1p2PI);
		}
		*AddKerfGCode = pNoKerfG;
		AddKerfGCode->R-=kerfvalue;
		AddKerfGCode->X0 += KerfX;
		AddKerfGCode->Y0 += KerfY;

		if(dir==G41)
		{
			KerfX = (kerfvalue)*cos(NextTangent+_1p2PI);
			KerfY = (kerfvalue)*sin(NextTangent+_1p2PI);	
		}
		else
		{
			KerfX = (kerfvalue)*cos(NextTangent-_1p2PI);
			KerfY = (kerfvalue)*sin(NextTangent-_1p2PI);
		}
		
		AddKerfGCode->X  += KerfX;
		AddKerfGCode->Y  += KerfY;
		AddKerfGCode->I = pNoKerfG.I;
		AddKerfGCode->J = pNoKerfG.J;
		
		
		
	}
	else if(pNoKerfG.Name==G03)
	{
		LastTangent = myatan2(pNoKerfG.Y0-pNoKerfG.J, pNoKerfG.X0-pNoKerfG.I)+_1p2PI;
		NextTangent = myatan2(pNoKerfG.Y-pNoKerfG.J, pNoKerfG.X-pNoKerfG.I)+_1p2PI;
		if(dir==G41)
		{
			KerfX = (kerfvalue)*cos(LastTangent+_1p2PI);
			KerfY = (kerfvalue)*sin(LastTangent+_1p2PI);
		}
		else
		{
			KerfX = (kerfvalue)*cos(LastTangent-_1p2PI);
			KerfY = (kerfvalue)*sin(LastTangent-_1p2PI);
		}
		*AddKerfGCode = pNoKerfG;
		AddKerfGCode->R+=kerfvalue;
		AddKerfGCode->X0 += KerfX;
		AddKerfGCode->Y0 += KerfY;
		if(dir==G41)
		{
			KerfX = (kerfvalue)*cos(NextTangent+_1p2PI);
			KerfY = (kerfvalue)*sin(NextTangent+_1p2PI);
		}
		else
		{
			KerfX = (kerfvalue)*cos(NextTangent-_1p2PI);
			KerfY = (kerfvalue)*sin(NextTangent-_1p2PI);
		}
		AddKerfGCode->X  += KerfX;
		AddKerfGCode->Y  += KerfY;
		AddKerfGCode->I = pNoKerfG.I;
		AddKerfGCode->J = pNoKerfG.J;
	}
	return 0;
}

/*
	前一条线的末端切线角度 beta，后一条线的初端切线角度alpha
	alpha-beta>0 左补偿为截取，右补偿为增加G03
	alpha-beta<0 左补偿为增加G02，右补偿为截取
	alpha-beta ~=0 则根据以原则截取或增加不直线段 
	kerfvalue < 0为左补偿，
	kerfvalue>0为右补偿
	pPreviousLine未补偿的前一句G 代码
	pNextLine 未补偿的后一句G 代码
	AddLine补偿的G代码, Name==G01(补偿的是小直线 ),G02(左补偿),G03(左补偿),M02(截取了，没有补偿)
	
*/
int Kerf::AddOrTrunc(GCodeARRAY_STRUCT &pPreviousLine,GCodeARRAY_STRUCT &pNextLine,GCodeARRAY_STRUCT *pAddLine,double kerfvalue,int dir)
{
	double AngleChange=0;
	//long InterTrunc;
	
	AngleChange = GetTangent(&pNextLine,0)-GetTangent(&pPreviousLine,1);
	
	while(AngleChange>(PI))
	{
		AngleChange-=2*PI;
	}
	while(AngleChange<-PI)
	{
		AngleChange+=2*PI;
	}
	if(IsLesser( fabs( fabs(AngleChange)-PI ),0.017))//1//1degree
	{
		if(dir == G41) //左补偿
			pAddLine->Name = G02;
		else
			pAddLine->Name = G03;	
	}
	else if(IsLesser(fabs(AngleChange),0.0055))///1 degree
	{
		pAddLine->Name = M00;
	}
	/*else if(IsLesser(fabs(AngleChange),0.035))
	{
		pAddLine->Name = G01;	
	}*/
	else if(dir == G41) //左补偿
	{
		if(AngleChange>0)
			pAddLine->Name = M02;
		else
		{	
			if(IsLesser(fabs(AngleChange),0.085))
				pAddLine->Name = G01;	
			else
				pAddLine->Name = G02;	
		}
	}
	else
	{
		if(AngleChange<0)
			pAddLine->Name = M02;
		else
		{
			if(IsLesser(fabs(AngleChange),0.017))
				pAddLine->Name = G01;	
			else
				pAddLine->Name = G03;	
		}
	}
  return 0;
}
int Kerf::Setupkerf(GCodeARRAY_STRUCT *pGcode, double *dx,double *dy,double kerfvlaue, int dir)
{
	double QieXianAngle1;
	if(pGcode->Name==G01 ||pGcode->Name==G00)
	{
		QieXianAngle1 = myatan2(pGcode->Y-pGcode->Y0,pGcode->X-pGcode->X0);
	}
	else if(pGcode->Name==G02 )
	{
		QieXianAngle1 = myatan2(pGcode->Y0-pGcode->J,pGcode->X0-pGcode->I)-_1p2PI;	
	}
	else if(pGcode->Name==G03 )
	{
		QieXianAngle1 = myatan2(pGcode->Y0-pGcode->J,pGcode->X0-pGcode->I)+_1p2PI;
	}
	if(dir==G41)
	{
		*dx = kerfvlaue*cos(QieXianAngle1+_1p2PI);
		*dy = kerfvlaue*sin(QieXianAngle1+_1p2PI);
	}
	else
	{
		*dx = kerfvlaue*cos(QieXianAngle1-_1p2PI);
		*dy = kerfvlaue*sin(QieXianAngle1-_1p2PI);
	}
  return 0;
}

int Kerf::Canclekerf(GCodeARRAY_STRUCT *pGcode, double *dx,double *dy,double kerfvlaue, int dir)
{
	double QieXianAngle1;
	if(pGcode->Name==G01||pGcode->Name==G00 )
	{
		QieXianAngle1 = myatan2(pGcode->Y-pGcode->Y0,pGcode->X-pGcode->X0);
	}
	else if(pGcode->Name==G02 )
	{
		QieXianAngle1 = myatan2(pGcode->Y-pGcode->J,pGcode->X-pGcode->I)-_1p2PI;	
	}
	else if(pGcode->Name==G03 )
	{
		QieXianAngle1 = myatan2(pGcode->Y-pGcode->J,pGcode->X-pGcode->I)+_1p2PI;
	}
	if(dir==G41)
	{
		*dx = kerfvlaue*cos(QieXianAngle1-_1p2PI);
		*dy = kerfvlaue*sin(QieXianAngle1-_1p2PI);
	}
	else
	{
		*dx = kerfvlaue*cos(QieXianAngle1+_1p2PI);
		*dy = kerfvlaue*sin(QieXianAngle1+_1p2PI);
	}
  return 0;
}

//#define  Assert(p) GUI_DispStringAt(p,col[30],row[10])
void Kerf::g2kerf(std::vector<GCodeARRAY_STRUCT> &DesKerFile,
                  std::vector<GCodeARRAY_STRUCT> &NoKerfFile) {

	enum kerfstatus
	{
		NOKERF,  //未建立 补偿
		G41KERF,  //遇到G41，但补偿尚没有建立状态 
		G42KERF,  //遇到G42,但补偿尚未建立状态
		JUSTSETG41KERF,//刚刚建立 G41补偿
		JUSTSETG42KERF,//刚刚建立 g42补偿
		SETTEDG41KERF, //已经建立起G41补偿
		SETTEDG42KERF, //已经建立起G42补偿
		G40KERF  //撤消补偿状态
	};
	std::vector<GCodeARRAY_STRUCT>::iterator GCodeArryPtrSrc = NoKerfFile.begin();
	std::vector<GCodeARRAY_STRUCT>::iterator GCodeArryPtrDes = DesKerFile.begin();
	GCodeARRAY_STRUCT GTemp1,GTemp2;
	char kerf_on=0;  // 1-建立割缝 
	double deltax,deltay,point2x,point2y;
	double x0,y0,x1,y1,cx,cy;
	int res,i;
	unsigned short LastGName,LastKerf=0,RowNum=0;
//	if(graphylimitxy.MaxRadius>50000)//半径大于50000，圆弧拆分
//	{
//		CircleCheFen(GfileFloatNoKerf);
//	}
	//UseKerf = 1;
/*	while(1)
    {
        if(GCodeArryPtrSrc->Name == G41 || GCodeArryPtrSrc->Name==G42)
            res = 1;
        if(GCodeArryPtrSrc->Name == G01 ||GCodeArryPtrSrc->Name == G02||GCodeArryPtrSrc->Name == G03)
        {
            if(res==1)
                UseKerf = 1;// 代码中有割缝补偿，否则没有割缝补偿
            break;
        }
        else if(GCodeArryPtrSrc->Name == M02)
        	break;
        GCodeArryPtrSrc++;
       // return ;
    } 
*/
  GCodeArryPtrSrc = GfileFloatNoKerf.begin();
 	
	if(((kerf <= EP) && (kerf >= -EP))/*||UseKerf==0*/)//if no kerf or zero kerf
	{
		while(GCodeArryPtrSrc->Name!=M02)
    {   
        GfileFloatKerf.push_back(*GCodeArryPtrSrc++);
    }
			
		GfileFloatKerf.push_back(*GCodeArryPtrSrc);
		return;
	}
	GCodeArryPtrSrc = GfileFloatNoKerf.begin();
  GfileFloatKerf.push_back(*GCodeArryPtrSrc);
  GCodeArryPtrDes = GfileFloatKerf.begin();
	while(GCodeArryPtrSrc->Name !=M02)
	{
		
		if(GCodeArryPtrSrc->Name==G41||GCodeArryPtrSrc->Name==G42 )
		{
			LastGName = GCodeArryPtrSrc->Name;
			if(kerf_on!=SETTEDG41KERF&&kerf_on!=SETTEDG42KERF)
			{
				i=0;
				do{
					
					GTemp1 = *(GCodeArryPtrSrc+i);
					i++;
					if(GTemp1.Name==M02||GTemp1.Name==G00||GTemp1.Name==G01  ||GTemp1.Name==G02 || GTemp1.Name==G03)
					break;
				} while(1);
				
				
	       if(GTemp1.Name!=M02)
	       {
	       	Setupkerf(&GTemp1, &deltax, &deltay, kerf, LastGName);
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *GCodeArryPtrSrc;
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G00;  //增加一段空走线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GCodeArryPtrDes->X0+deltax;
					GCodeArryPtrDes->Y = GCodeArryPtrDes->Y0+deltay;
        }
				if(LastGName==G41)
					kerf_on = JUSTSETG41KERF;
				else if(LastGName==G42)
					kerf_on = JUSTSETG42KERF;
			}
		}
		else if(GCodeArryPtrSrc->Name==G40)  //取消补偿
		{
			if(kerf_on==SETTEDG41KERF||kerf_on==SETTEDG42KERF)
			{
				if(kerf_on == SETTEDG41KERF)
					LastGName = G41;
				else
					LastGName = G42;
				i=0;
				do {
					GTemp1 = *(GCodeArryPtrDes-i);
					i++;
				} while (GTemp1.ShowLine>1 && GTemp1.Name!=G01  && GTemp1.Name!=G02 && GTemp1.Name!=G03);

				Canclekerf(&GTemp1, &deltax, &deltay, kerf, LastGName);
			  //GCodeArryPtrDes++;
			  //*GCodeArryPtrDes = *GCodeArryPtrSrc;
        GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
        GCodeArryPtrDes = GfileFloatKerf.end() - 1;
				GCodeArryPtrDes->Name = G00;  //增加一段空走线
				GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
				GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
				GCodeArryPtrDes->X = GCodeArryPtrDes->X0+deltax;
				GCodeArryPtrDes->Y = GCodeArryPtrDes->Y0+deltay;
				kerf_on = NOKERF;
				LastKerf = NOKERF;
			}
		}
		else if(GCodeArryPtrSrc->Name==G01&&(fabs(GCodeArryPtrSrc->X-GCodeArryPtrSrc->X0)>0.01||fabs(GCodeArryPtrSrc->Y-GCodeArryPtrSrc->Y0)>0.01))  // case G01 level 1
		{
			
			if(kerf_on==JUSTSETG41KERF ||kerf_on==JUSTSETG42KERF)
			{
				if(kerf_on==JUSTSETG42KERF)
				{
					LastGName = G42;
					kerf_on = SETTEDG42KERF;
				}
				else
				{
					LastGName = G41;
					kerf_on = SETTEDG41KERF;
				}
				
				GetAddKerfGCode(*GCodeArryPtrSrc, &GTemp1, kerf,LastGName);
				//GCodeArryPtrDes++;
				//*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
        GfileFloatKerf.push_back(*GCodeArryPtrDes); // Add by ZhanShi
	      GCodeArryPtrDes = GfileFloatKerf.end() - 1;
				GCodeArryPtrDes->Name = G01;
				GCodeArryPtrDes->X0 = GTemp1.X0;
				GCodeArryPtrDes->Y0 = GTemp1.Y0;
				GCodeArryPtrDes->X = GTemp1.X;
				GCodeArryPtrDes->Y = GTemp1.Y;
			}
			else if(kerf_on == SETTEDG41KERF||kerf_on==SETTEDG42KERF)  //仅仅确定当前行的起点坐标，终点坐标采用补偿后的末点坐标
			{
				if(kerf_on==SETTEDG42KERF)
				{
					AddOrTrunc(*(GCodeArryPtrSrc-1), *GCodeArryPtrSrc, &GTemp1, kerf, G42);
					GetAddKerfGCode(*GCodeArryPtrSrc, &GTemp2, kerf,G42);
				}
				else
				{
					AddOrTrunc(*(GCodeArryPtrSrc-1), *GCodeArryPtrSrc, &GTemp1, kerf, G41);
					GetAddKerfGCode(*GCodeArryPtrSrc, &GTemp2, kerf,G41);
				}
				if(GTemp1.Name==M02)  //前面一行和当前行需要截取
				{
					if((GCodeArryPtrSrc-1)->Name==G01)  //case G01 level 2
					{
						res = TwoLineIsIntersect(GCodeArryPtrDes->X0,GCodeArryPtrDes->Y0, GCodeArryPtrDes->X,GCodeArryPtrDes->Y, GTemp2.X0, GTemp2.Y0,GTemp2.X, GTemp2.Y, &deltax,&deltay);
						if(res==1)
						{
							GCodeArryPtrDes->X = deltax;
							GCodeArryPtrDes->Y = deltay;
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
              GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;
							GCodeArryPtrDes->X0 = deltax;
							GCodeArryPtrDes->Y0 = deltay;
							GCodeArryPtrDes->X = GTemp2.X;
							GCodeArryPtrDes->Y = GTemp2.Y;
						
						}
						else if(res==0)
						{
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
              GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X;
							GCodeArryPtrDes->Y = GTemp2.Y;
							//Assert("line line res=0");
							//_sleep(20);
						}
						else
						{
						//没有交点，则寻求与前一行的交点
							/*i=1;
							do{
								*pTemp = *(GCodeArryPtrDes-i);
								i++;
								}while(pTemp->ShowLine>1&&pTemp->Name!=G01&&pTemp->Name!=G02&&pTemp->Name!=G03);
							if(pTemp->Name==G01)
							{
								res = TwoLineIsIntersect(GTemp2.X0, GTemp2.Y0, GTemp2.X, GTemp2.Y, pTemp->X0, pTemp->Y0, pTemp->X, pTemp->Y, &point2x, &point2y);
								if(res==1)
								{
									pTemp->X = deltax;
									pTemp->Y = deltay;
									//GCodeArryPtrDes++;
									GCodeArryPtrDes = pTemp+1;
									*GCodeArryPtrDes = *(GCodeArryPtrSrc);
									GCodeArryPtrDes->Name = G01;
									GCodeArryPtrDes->X0 = deltax;
									GCodeArryPtrDes->Y0 = deltay;
									GCodeArryPtrDes->X = GTemp2.X;
									GCodeArryPtrDes->Y = GTemp2.Y;
								}
								else if(res==0)
								{
									GCodeArryPtrDes++;
									*GCodeArryPtrDes = *(GCodeArryPtrSrc);
									GCodeArryPtrDes->Name = G01;
									GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
									GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
									GCodeArryPtrDes->X = GTemp2.X;
									GCodeArryPtrDes->Y = GTemp2.Y;
									Assert("ErrorG01/G01/G01");
								}
								
							}
							else if(pTemp->Name==G02 || pTemp->Name==G03)
							{
								res = LineCircleIntersect(GTemp2.X0, GTemp2.Y0, GTemp2.X, GTemp2.Y, pTemp->I, pTemp->J, GetRadius(pTemp), &deltax, &deltay, &point2x,&point2y);
								if(res==0) //直线与圆没有交点，则插入小直线
								{
									Assert("G01/G01/G0203");
								}
								else if(res==1)

								{
									pTemp->X = deltax;
									pTemp->Y = deltay;
									GCodeArryPtrDes=pTemp+1;
									*GCodeArryPtrDes = *(GCodeArryPtrSrc);
									GCodeArryPtrDes->Name = G01;
									GCodeArryPtrDes->X0 = deltax;
									GCodeArryPtrDes->Y0 = deltay;
									GCodeArryPtrDes->X = GTemp2.X;
									GCodeArryPtrDes->Y = GTemp2.Y;
								}
								else if(res==2)
								{
									if(point_in_line(point2x,point2y, GTemp2.X0,GTemp2.Y0,GTemp2.X,GTemp2.Y))
									//if(dist(point2x,point2y, GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0)<dist(deltax,deltay, GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0))
									{
										pTemp->X = point2x;
										pTemp->Y = point2y;
										GCodeArryPtrDes=pTemp+1;
										*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
										GCodeArryPtrDes->Name = G01;
										GCodeArryPtrDes->X0 = point2x;
										GCodeArryPtrDes->Y0 = point2y;
										GCodeArryPtrDes->X = GTemp2.X;
										GCodeArryPtrDes->Y = GTemp2.Y;
									}
									else
									{
										pTemp->X = deltax;
										pTemp->Y = deltay;
										GCodeArryPtrDes=pTemp+1;
										*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
										GCodeArryPtrDes->Name = G01;
										GCodeArryPtrDes->X0 = deltax;
										GCodeArryPtrDes->Y0 = deltay;
										GCodeArryPtrDes->X = GTemp2.X;
										GCodeArryPtrDes->Y = GTemp2.Y;
									}
								}
							}*/
							
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
              GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X;
							GCodeArryPtrDes->Y = GTemp2.Y;
							//Assert("line line no joint");
							//_sleep(20);
							//assert(0);

						}
					}
					else if((GCodeArryPtrSrc-1)->Name==G02||(GCodeArryPtrSrc-1)->Name==G03) //case G01 level 2,上行G02或G03，当前行G01
					{
						res = LineCircleIntersect(GTemp2.X0, GTemp2.Y0, GTemp2.X, GTemp2.Y, GCodeArryPtrDes->I, GCodeArryPtrDes->J, GetRadius(*GCodeArryPtrDes), &deltax, &deltay, &point2x,&point2y);
						if(res==0) //直线与圆没有交点，则插入小直线
						{
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
              GfileFloatKerf.push_back(*GCodeArryPtrDes); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X0;
							GCodeArryPtrDes->Y = GTemp2.Y0;
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
              GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X;
							GCodeArryPtrDes->Y = GTemp2.Y;
							//Assert("C/L no joint");
							//_sleep(20);
							//assert(0);

						}
						else if(res==1)
						{
							if(point_in_line(deltax,deltay,GTemp2.X0, GTemp2.Y0,GTemp2.X, GTemp2.Y))
							{
								GCodeArryPtrDes->X = deltax;
								GCodeArryPtrDes->Y = deltay;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = G01;
								GCodeArryPtrDes->X0 = deltax;
								GCodeArryPtrDes->Y0 = deltay;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;

							}
							else
							{
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
                GfileFloatKerf.push_back(*GCodeArryPtrDes); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = G01;
								GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
								GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
								GCodeArryPtrDes->X = GTemp2.X0;
								GCodeArryPtrDes->Y = GTemp2.Y0;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = G01;
								GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
								GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								//Assert("C/L apart");
								//_sleep(20);
								//assert(0);

							}
								
						}
						else
						{
							if(dist(point2x,point2y, GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0)<dist(deltax,deltay, GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0))
							{
								GCodeArryPtrDes->X = point2x;
								GCodeArryPtrDes->Y = point2y;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = G01;
								GCodeArryPtrDes->X0 = point2x;
								GCodeArryPtrDes->Y0 = point2y;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
							}
							else
							{
								/*if((GCodeArryPtrSrc-1)->Name==G02)
								{
									x0 = (GCodeArryPtrSrc-1)->X;
									y0 = (GCodeArryPtrSrc-1)->Y;
									x1 = (GCodeArryPtrSrc-1)->X0;
									y1 = (GCodeArryPtrSrc-1)->Y0;
									cx = (GCodeArryPtrSrc-1)->I;
									cy = (GCodeArryPtrSrc-1)->J;
								}
								else
								{
									x0 = (GCodeArryPtrSrc-1)->X0;
									y0 = (GCodeArryPtrSrc-1)->Y0;
									x1 = (GCodeArryPtrSrc-1)->X;
									y1 = (GCodeArryPtrSrc-1)->Y;
									cx = (GCodeArryPtrSrc-1)->I;
									cy = (GCodeArryPtrSrc-1)->J;
								} 2014.6.16 */
								if((GCodeArryPtrDes)->Name==G02) 
								{
									x0 = (GCodeArryPtrDes)->X;
									y0 = (GCodeArryPtrDes)->Y;
									x1 = (GCodeArryPtrDes)->X0;
									y1 = (GCodeArryPtrDes)->Y0;
									cx = (GCodeArryPtrDes)->I;
									cy = (GCodeArryPtrDes)->J;
								}
								else
								{
									x0 = (GCodeArryPtrDes)->X0;
									y0 = (GCodeArryPtrDes)->Y0;
									x1 = (GCodeArryPtrDes)->X;
									y1 = (GCodeArryPtrDes)->Y;
									cx = (GCodeArryPtrDes)->I;
									cy = (GCodeArryPtrDes)->J;
								}
								
								//?????此处先简单处理，截取没有的圆弧忽略，起点从上一行的末点
								if(if_point_in_arc(x0,y0,x1,y1,cx,cy,deltax,deltay))
								{
									GCodeArryPtrDes->X = deltax;
									GCodeArryPtrDes->Y = deltay;
									//GCodeArryPtrDes++;
									//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                  GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                  GCodeArryPtrDes = GfileFloatKerf.end() - 1;
									GCodeArryPtrDes->Name = G01;
									GCodeArryPtrDes->X0 = deltax;
									GCodeArryPtrDes->Y0 = deltay;
									GCodeArryPtrDes->X = GTemp2.X;
									GCodeArryPtrDes->Y = GTemp2.Y;
	
								}
								else
								{
									GCodeArryPtrDes->X = point2x;
									GCodeArryPtrDes->Y = point2y;
									//GCodeArryPtrDes++;
									//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                  GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                  GCodeArryPtrDes = GfileFloatKerf.end() - 1;
									GCodeArryPtrDes->Name = G01;
									GCodeArryPtrDes->X0 = point2x;
									GCodeArryPtrDes->Y0 = point2y;
									GCodeArryPtrDes->X = GTemp2.X;
									GCodeArryPtrDes->Y = GTemp2.Y;
									
									/*
									GCodeArryPtrDes = *(GCodeArryPtrSrc);
									GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
									GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
									GCodeArryPtrDes->X = GTemp2.X;
									GCodeArryPtrDes->Y = GTemp2.Y;
									GCodeArryPtrDes->Name = G01;*/
								}
							}
						}
					}
				}
				else if(GTemp1.Name==G01) //前面一行和当前行需要插入一条小直线 
				{
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G01;  //插入小直线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X0;
					GCodeArryPtrDes->Y = GTemp2.Y0;

          //确定下一条线的起点								
					//GCodeArryPtrDes++;    
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G01;  //插入小直线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
					
				}
				else if(GTemp1.Name==G02)
				{
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					if(sqrt(pow((GCodeArryPtrDes-1)->X-GTemp2.X0,2)+pow((GCodeArryPtrDes-1)->Y-GTemp2.Y0,2))<0.2)
					GCodeArryPtrDes->Name = G01;
					else
					GCodeArryPtrDes->Name = G02;  //插入小圆弧
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X0;
					GCodeArryPtrDes->Y = GTemp2.Y0;
					GCodeArryPtrDes->I =GCodeArryPtrSrc->X0;
					GCodeArryPtrDes->J =GCodeArryPtrSrc->Y0;

					//确定下一条线的起点			
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G01;  //插入小直线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
					

				}
				else if(GTemp1.Name==G03)
				{
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					if(sqrt(pow((GCodeArryPtrDes-1)->X-GTemp2.X0,2)+pow((GCodeArryPtrDes-1)->Y-GTemp2.Y0,2))<0.2)
					GCodeArryPtrDes->Name = G01;
					else
					GCodeArryPtrDes->Name = G03;  //插入小圆弧
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X0;
					GCodeArryPtrDes->Y = GTemp2.Y0;
					GCodeArryPtrDes->I =GCodeArryPtrSrc->X0;
					GCodeArryPtrDes->J =GCodeArryPtrSrc->Y0;
	
								
					//确定下一条线的起点			
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G01;  //插入小直线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
		
								
				}
				else if(GTemp1.Name==M00) //前面一行和当前行平等相接，不需要插入或截取
				{
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
					
								
				}

			}
			else
			{
				//GCodeArryPtrDes++;
				//*GCodeArryPtrDes = *GCodeArryPtrSrc;
        GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
        GCodeArryPtrDes = GfileFloatKerf.end() - 1;
				GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
				GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
				GCodeArryPtrDes->X = (GCodeArryPtrDes)->X0+GCodeArryPtrSrc->X - GCodeArryPtrSrc->X0;
				GCodeArryPtrDes->Y = (GCodeArryPtrDes)->Y0+GCodeArryPtrSrc->Y - GCodeArryPtrSrc->Y0;
				
			}
		}
		else if( GCodeArryPtrSrc->Name==G02 || GCodeArryPtrSrc->Name==G03)  //case G02/G03 level 1
		{
			/*if(kerf_on == G40KERF)  //防止有的代码只出现一次G41或G42，中间多次穿孔，最后才一个G40
			{
			//此段代码和if(GCodeArryPtrSrc->Name==G41||GCodeArryPtrSrc->Name==G42)段差不多
			       
				
				if(LastKerf==SETTEDG41KERF|LastKerf==SETTEDG42KERF)
				{
					i=0;
					GTemp1 =  *GCodeArryPtrSrc;
				  Setupkerf(&GTemp1, &deltax, &deltay, kerf, GCodeArryPtrSrc->Name);
					GCodeArryPtrDes++;
					*GCodeArryPtrDes = *GCodeArryPtrSrc;
					GCodeArryPtrDes->Name = G00;  //增加一段空走线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GCodeArryPtrDes->X0+deltax;
					GCodeArryPtrDes->Y = GCodeArryPtrDes->Y0+deltay;
					if(LastKerf==SETTEDG41KERF)
						kerf_on = JUSTSETG41KERF;
					else
						kerf_on = JUSTSETG42KERF;
					LastGName  =kerf_on;
				}
			}*/
			
			if(kerf_on==JUSTSETG41KERF ||kerf_on==JUSTSETG42KERF)
			{
				if(kerf_on==JUSTSETG42KERF)
				{
					LastGName = G42;
					kerf_on = SETTEDG42KERF;
				}
				else
				{
					LastGName = G41;
					kerf_on = SETTEDG41KERF;
				}
				
				GetAddKerfGCode(*GCodeArryPtrSrc, &GTemp1, kerf,LastGName);
				//GCodeArryPtrDes++;
				//*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
        GfileFloatKerf.push_back(*GCodeArryPtrDes); // Add by ZhanShi
        GCodeArryPtrDes = GfileFloatKerf.end() - 1;
				GCodeArryPtrDes->Name = GCodeArryPtrSrc->Name;
				GCodeArryPtrDes->X0 = GTemp1.X0;
				GCodeArryPtrDes->Y0 = GTemp1.Y0;
				GCodeArryPtrDes->X = GTemp1.X;
				GCodeArryPtrDes->Y = GTemp1.Y;
				GCodeArryPtrDes->I = GTemp1.I;
				GCodeArryPtrDes->J = GTemp1.J;

								
			}
			
			else if(kerf_on == SETTEDG41KERF||kerf_on==SETTEDG42KERF)  //仅仅确定当前行的起点坐标，终点坐标采用补偿后的末点坐标
			{
				LastGName = GCodeArryPtrSrc->Name;
				if(kerf_on==SETTEDG42KERF)
				{
					AddOrTrunc(*(GCodeArryPtrSrc-1), *GCodeArryPtrSrc, &GTemp1, kerf, G42);
					GetAddKerfGCode(*GCodeArryPtrSrc, &GTemp2, kerf,G42);
				}
				else
				{
					AddOrTrunc(*(GCodeArryPtrSrc-1), *GCodeArryPtrSrc, &GTemp1, kerf, G41);
					GetAddKerfGCode(*GCodeArryPtrSrc, &GTemp2, kerf,G41);
				}
				if(GTemp1.Name==M02)  //前面一行和当前行需要截取
        {
					if((GCodeArryPtrSrc-1)->Name==G01)
					{
						res = LineCircleIntersect(GCodeArryPtrDes->X0 ,GCodeArryPtrDes->Y0 ,GCodeArryPtrDes->X,GCodeArryPtrDes->Y, GTemp2.I, GTemp2.J, GetRadius(GTemp2), &deltax, &deltay, &point2x,&point2y);
						if(res==0)
						{
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
              GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;  //没有相交则插入小直线
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X0;
							GCodeArryPtrDes->Y = GTemp2.Y0;

							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = GTemp2;  //当前行
              GfileFloatKerf.push_back(GTemp2); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X;
							GCodeArryPtrDes->Y = GTemp2.Y;
							GCodeArryPtrDes->I = GTemp2.I;
							GCodeArryPtrDes->J = GTemp2.J;
							
								
							//Assert("res=0,line and circle") ;
              //_sleep(20);
							//assert(0);

								
						}
						else if(res==1)
						{
							if(point_in_line(deltax,deltay,(GCodeArryPtrSrc-1)->X0, (GCodeArryPtrSrc-1)->Y0,(GCodeArryPtrSrc-1)->X, (GCodeArryPtrSrc-1)->Y))
							{
								GCodeArryPtrDes->X = deltax;
								GCodeArryPtrDes->Y = deltay;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = deltax;
								GCodeArryPtrDes->Y0 = deltay;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

								
								
							}
							else
							{
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = G01;  //没有相交则插入小直线
								GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
								GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
								GCodeArryPtrDes->X = GTemp2.X0;
								GCodeArryPtrDes->Y = GTemp2.Y0;
	
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = GTemp2;  //当前行
                GfileFloatKerf.push_back(GTemp2); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
								GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;
								
								//Assert("res=1,L/C apart") ;
								//_sleep(20);
								//assert(0);

							}
								
						}
						else
						{
							if(dist(point2x,point2y,GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0)<dist(deltax,deltay,GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0))
							{
								GCodeArryPtrDes->X = point2x;
								GCodeArryPtrDes->Y = point2y;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = point2x;
								GCodeArryPtrDes->Y0 = point2y;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

							}
							else
							{
								GCodeArryPtrDes->X = deltax;
								GCodeArryPtrDes->Y = deltay;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = deltax;
								GCodeArryPtrDes->Y0 = deltay;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

							}
							
						}
					}
					else if((GCodeArryPtrSrc-1)->Name==G02||(GCodeArryPtrSrc-1)->Name==G03)
					{
						res = CircleIntersector(GCodeArryPtrDes->I, GCodeArryPtrDes->J, GetRadius(*GCodeArryPtrDes),GTemp2.I, GTemp2.J,GetRadius(GTemp2), &deltax, &deltay, &point2x,&point2y);
						if(res==1 || res == 5) //相离 或 相含
						{
							//增加小直线
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrDes-1);
              GfileFloatKerf.push_back(*GCodeArryPtrDes); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = G01;  //插入小圆弧
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X0;
							GCodeArryPtrDes->Y = GTemp2.Y0;
	
							//确定下一条线的起点	
							//GCodeArryPtrDes++;
							//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
              GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
              GCodeArryPtrDes = GfileFloatKerf.end() - 1;
							GCodeArryPtrDes->Name = LastGName;  //抄过当前圆弧
							GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
							GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
							GCodeArryPtrDes->X = GTemp2.X;
							GCodeArryPtrDes->Y = GTemp2.Y;
							GCodeArryPtrDes->I = GTemp2.I;
							GCodeArryPtrDes->J = GTemp2.J;
							//Assert("C/C apart");
						  //_sleep(20);
							//assert(0);

						}
						else if(res==2||res==3)
						{
							//离原点近的交点为正确的交点
							if(dist(GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0,deltax,deltay)<dist(GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0,point2x,point2y))
							{
								GCodeArryPtrDes->X = deltax;
								GCodeArryPtrDes->Y = deltay;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = deltax;
								GCodeArryPtrDes->Y0 = deltay;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

							}
							else
							{
								GCodeArryPtrDes->X = point2x;
								GCodeArryPtrDes->Y = point2y;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = deltax;
								GCodeArryPtrDes->Y0 = deltay;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

							}
								
						}
						else
						{
							//离原点近的交点为正确的交点
							if(dist(GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0,deltax,deltay)<dist(GCodeArryPtrSrc->X0,GCodeArryPtrSrc->Y0,point2x,point2y))
							{
								GCodeArryPtrDes->X = deltax;
								GCodeArryPtrDes->Y = deltay;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = deltax;
								GCodeArryPtrDes->Y0 = deltay;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

							}
							else
							{
								GCodeArryPtrDes->X = point2x;
								GCodeArryPtrDes->Y = point2y;
								//GCodeArryPtrDes++;
								//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
                GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
                GCodeArryPtrDes = GfileFloatKerf.end() - 1;
								GCodeArryPtrDes->Name = LastGName;
								GCodeArryPtrDes->X0 = point2x;
								GCodeArryPtrDes->Y0 = point2y;
								GCodeArryPtrDes->X = GTemp2.X;
								GCodeArryPtrDes->Y = GTemp2.Y;
								GCodeArryPtrDes->I = GTemp2.I;
								GCodeArryPtrDes->J = GTemp2.J;

							}
							
						}
					}
				}
				else if(GTemp1.Name==G01) //前面一行和当前行需要插入一条小直线 
        {
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G01;  //插入小直线
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X0;
					GCodeArryPtrDes->Y = GTemp2.Y0;

					//确定下一条线的起点			
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = LastGName;  //插入当前G02 或G03
					GCodeArryPtrDes->X0 = GTemp2.X0;
					GCodeArryPtrDes->Y0 = GTemp2.Y0;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
					GCodeArryPtrDes->I = GTemp2.I;
					GCodeArryPtrDes->J = GTemp2.J;

					
				}
				else if(GTemp1.Name==G02)
        {
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G02;  //插入小圆弧
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X0;
					GCodeArryPtrDes->Y = GTemp2.Y0;
					GCodeArryPtrDes->I =GCodeArryPtrSrc->X0;
					GCodeArryPtrDes->J =GCodeArryPtrSrc->Y0;
	
								
					//确定下一条线的起点
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = LastGName;  //插入当前G02或G03
					GCodeArryPtrDes->X0 = GTemp2.X0;
					GCodeArryPtrDes->Y0 = GTemp2.Y0;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
					GCodeArryPtrDes->I = GTemp2.I;
					GCodeArryPtrDes->J = GTemp2.J;

								
				}
				else if(GTemp1.Name==G03)
        {
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = G03;  //插入小圆弧
					GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
					GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
					GCodeArryPtrDes->X = GTemp2.X0;
					GCodeArryPtrDes->Y = GTemp2.Y0;
					GCodeArryPtrDes->I =GCodeArryPtrSrc->X0;
					GCodeArryPtrDes->J =GCodeArryPtrSrc->Y0;

								
					//确定下一条线的起点
					//GCodeArryPtrDes++;
					//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
          GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
					GCodeArryPtrDes->Name = LastGName;  //插入当前G02或G03
					GCodeArryPtrDes->X0 = GTemp2.X0;
					GCodeArryPtrDes->Y0 = GTemp2.Y0;
					GCodeArryPtrDes->X = GTemp2.X;
					GCodeArryPtrDes->Y = GTemp2.Y;
					GCodeArryPtrDes->I = GTemp2.I;
					GCodeArryPtrDes->J = GTemp2.J;

								
				}
				else if(GTemp1.Name==M00) //前面一行和当前行平等相接，不需要插入或截取
        {
					if (fabs((GCodeArryPtrDes)->X-GTemp2.X0)>0.1||
					    fabs((GCodeArryPtrDes)->Y-GTemp2.Y0)>0.1)
					{
						//需要插入直线
            GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
            GCodeArryPtrDes = GfileFloatKerf.end() - 1;
						GCodeArryPtrDes->Name = G01;  
						GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
						GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
						GCodeArryPtrDes->X = GTemp2.X0;
						GCodeArryPtrDes->Y = GTemp2.Y0;
						
						//再把当前圆弧加入进来
            //GCodeArryPtrDes++;
						//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
            GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
	          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
						GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
						GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
						GCodeArryPtrDes->X = GTemp2.X;
						GCodeArryPtrDes->Y = GTemp2.Y;
						GCodeArryPtrDes->I = GTemp2.I;
						GCodeArryPtrDes->J = GTemp2.J;

					
					}
					else
					{
            //GCodeArryPtrDes++;
						//*GCodeArryPtrDes = *(GCodeArryPtrSrc);
            GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
	          GCodeArryPtrDes = GfileFloatKerf.end() - 1;
						GCodeArryPtrDes->X0 = GTemp2.X0;//不取上一行的末点坐标，是因为有时起点坐标偏一点就会导致一个整圆变成小圆弧
						GCodeArryPtrDes->Y0 = GTemp2.Y0;
						GCodeArryPtrDes->X = GTemp2.X;
						GCodeArryPtrDes->Y = GTemp2.Y;
						GCodeArryPtrDes->I = GTemp2.I;
						GCodeArryPtrDes->J = GTemp2.J;
					}
				}

			}
			else
			{
				//GCodeArryPtrDes++;
				//*GCodeArryPtrDes = *GCodeArryPtrSrc;
        GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
	      GCodeArryPtrDes = GfileFloatKerf.end() - 1;
				GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
				GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
				GCodeArryPtrDes->X = GCodeArryPtrDes->X0+(GCodeArryPtrSrc->X-GCodeArryPtrSrc->X0);
				GCodeArryPtrDes->Y = GCodeArryPtrDes->Y0+(GCodeArryPtrSrc->Y-GCodeArryPtrSrc->Y0);
				GCodeArryPtrDes->I = GCodeArryPtrSrc->I;
				GCodeArryPtrDes->J = GCodeArryPtrSrc->J;

			}
		}
		
		else
		{
			if(GCodeArryPtrSrc == GfileFloatNoKerf.begin())
			{
				*GCodeArryPtrDes = *GCodeArryPtrSrc;
			}
			else
			{
				//GCodeArryPtrDes++;
				//*GCodeArryPtrDes = *GCodeArryPtrSrc;
        GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
        GCodeArryPtrDes = GfileFloatKerf.end() - 1;
				GCodeArryPtrDes->X0 = (GCodeArryPtrDes-1)->X;
				GCodeArryPtrDes->Y0 = (GCodeArryPtrDes-1)->Y;
				GCodeArryPtrDes->X = GCodeArryPtrDes->X0+(GCodeArryPtrSrc->X-GCodeArryPtrSrc->X0);
				GCodeArryPtrDes->Y = GCodeArryPtrDes->Y0+(GCodeArryPtrSrc->Y-GCodeArryPtrSrc->Y0);
				GCodeArryPtrDes->I = GCodeArryPtrSrc->I;
				GCodeArryPtrDes->J = GCodeArryPtrSrc->J;
			}

		}
    if((GCodeArryPtrDes->Name>=G00&&GCodeArryPtrDes->Name<=G03)||(GCodeArryPtrDes->Name>=G26&&GCodeArryPtrDes->Name<=G28))
    {
        RowNum++;
        GCodeArryPtrDes->PierceHoleNum = RowNum;
    }
		GCodeArryPtrSrc++;
	}
  GfileFloatKerf.push_back(*GCodeArryPtrSrc); // Add by ZhanShi
	GCodeArryPtrDes = GfileFloatKerf.end() - 1;
	GCodeArryPtrDes->Name = M02;
}

bool Kerf::ReadGCode(const std::string &file_name, std::vector<std::string> &code_lines) {
  std::ifstream in_file(file_name.c_str());
  if (!in_file) {
    return false;
  }
  std::string line_text = "";
  while (getline(in_file, line_text)) {
    code_lines.push_back(line_text);
  }

  in_file.close();
  in_file.clear();

  return true;
}

// assume that is an absolute g code.
void Kerf::GCodeParse(const std::vector<std::string> &code_lines) {
  GfileFloatNoKerf.clear();

  double start_X = 0.;
  double start_Y = 0.;

  for (size_t i = 0; i < code_lines.size(); i++) {
    GCodeARRAY_STRUCT code_array;
    code_array.X0 = code_array.X = start_X;
    code_array.Y0 = code_array.Y = start_Y;

    std::string code_type = code_lines[i].substr(0, 3);
    if (code_type.compare("G00") == 0) {
      code_array.Name = G00;
      code_array.X = is_absolute ? GetCodeValue(code_lines[i], "X") : start_X + GetCodeValue(code_lines[i], "X");
      code_array.Y = is_absolute ? GetCodeValue(code_lines[i], "Y") : start_Y + GetCodeValue(code_lines[i], "Y");

      code_array.X0 = start_X;
      code_array.Y0 = start_Y;
      start_X = code_array.X;
      start_Y = code_array.Y;

      code_array.Length = sqrt(pow(code_array.X - code_array.X0, 2)
          + pow(code_array.Y - code_array.Y0, 2));

    } else if (code_type.compare("G01") == 0) {
      code_array.Name = G01;
      code_array.X = is_absolute ? GetCodeValue(code_lines[i], "X") : start_X + GetCodeValue(code_lines[i], "X");
      code_array.Y = is_absolute ? GetCodeValue(code_lines[i], "Y") : start_Y + GetCodeValue(code_lines[i], "Y");
      code_array.F = GetCodeValue(code_lines[i], "F");

      code_array.X0 = start_X;
      code_array.Y0 = start_Y;
      start_X = code_array.X;
      start_Y = code_array.Y;

      code_array.Length = sqrt(pow(code_array.X - code_array.X0, 2)
          + pow(code_array.Y - code_array.Y0, 2));

    } else if (code_type.compare("G02") == 0) {
      code_array.Name = G02;
      code_array.X = is_absolute ? GetCodeValue(code_lines[i], "X") : start_X + GetCodeValue(code_lines[i], "X");
      code_array.Y = is_absolute ? GetCodeValue(code_lines[i], "Y") : start_Y + GetCodeValue(code_lines[i], "Y");
      code_array.I = is_absolute ? GetCodeValue(code_lines[i], "I") : start_X + GetCodeValue(code_lines[i], "I");
      code_array.J = is_absolute ? GetCodeValue(code_lines[i], "J") : start_Y + GetCodeValue(code_lines[i], "J");
      code_array.F = GetCodeValue(code_lines[i], "F");

      code_array.X0 = start_X;
      code_array.Y0 = start_Y;
      start_X = code_array.X;
      start_Y = code_array.Y;

      code_array.R = GetRadius(code_array);
      code_array.Length = 2 * code_array.R;
    } else if (code_type.compare("G03") == 0) {
      code_array.Name = G03;
      code_array.X = is_absolute ? GetCodeValue(code_lines[i], "X") : start_X + GetCodeValue(code_lines[i], "X");
      code_array.Y = is_absolute ? GetCodeValue(code_lines[i], "Y") : start_Y + GetCodeValue(code_lines[i], "Y");
      code_array.I = is_absolute ? GetCodeValue(code_lines[i], "I") : start_X + GetCodeValue(code_lines[i], "I");
      code_array.J = is_absolute ? GetCodeValue(code_lines[i], "J") : start_Y + GetCodeValue(code_lines[i], "J");
      code_array.F = GetCodeValue(code_lines[i], "F");

      code_array.X0 = start_X;
      code_array.Y0 = start_Y;
      start_X = code_array.X;
      start_Y = code_array.Y;

      code_array.R = GetRadius(code_array);
      code_array.Length = 2 * code_array.R;
    } else if (code_type.compare("G41") == 0) {
      code_array.Name = G41;
    } else if (code_type.compare("G42") == 0) {
      code_array.Name = G42;
    } else if (code_type.compare("G40") == 0) {
      code_array.Name = G40;
    } else if (code_type.compare("G04") == 0) {
      code_array.Name = G04;
      code_array.Delay = GetCodeValue(code_lines[i], "F");
    } else if (code_type.compare("G20") == 0) {
      code_array.Name = G20;
    } else if (code_type.compare("G21") == 0) {
      code_array.Name = G21;
    } else if (code_type.compare("G90") == 0) {
      code_array.Name = G90;
      is_absolute = true;
    } else if (code_type.compare("G91") == 0) {
      code_array.Name = G91;
      is_absolute = false;
    } else if (code_type.compare("G99") == 0) {
      code_array.Name = G99;
      code_array.ScaleFactor = GetCodeValue(code_lines[i], "X");
      code_array.RotateAngle = GetCodeValue(code_lines[i], "Y");
      code_array.HorizonMirror = GetCodeValue(code_lines[i], "I");
      code_array.VerticalMirror = GetCodeValue(code_lines[i], "J");
    } else if (code_type.compare("M00") == 0) {
      code_array.Name = M00;
    } else if (code_type.compare("M02") == 0) {
      code_array.Name = M02;
    } else if (code_type.compare("M07") == 0) {
      code_array.Name = M07;
    } else if (code_type.compare("M08") == 0) {
      code_array.Name = M08;
    } else if (code_type.compare("M09") == 0) {
      code_array.Name = M09;
    } else if (code_type.compare("M10") == 0) {
      code_array.Name = M10;
    } else if (code_type.compare("M11") == 0) {
      code_array.Name = M11;
    } else if (code_type.compare("M12") == 0) {
      code_array.Name = M12;
    }

    GfileFloatNoKerf.push_back(code_array);
  }

  GfileFloatKerf.clear();
  g2kerf(GfileFloatKerf, GfileFloatNoKerf);
}

double Kerf::GetCodeValue(const std::string &code_line,
                          const std::string &match) {

  std::size_t pos = code_line.find(match);
  if (pos == std::string::npos) {
    return 0.;
  }
  std::string value = code_line.substr(pos + 1);
  char *endptr;
  return strtod(value.c_str(), &endptr);
}

bool Kerf::WriteGCode(const std::string &file_name,
                      const std::vector<std::string> &contents) {

  std::ofstream out_file(file_name.c_str());
  if (!out_file) {
    return false;
  }
  for (size_t i = 0; i < contents.size(); i++) {
    out_file << contents[i] << "\n";
  }
  out_file.close();
  out_file.clear();
  return true;
}

void Kerf::GenerateKerfGCode(std::vector<std::string> &code_lines) {
  code_lines.clear();
  for (size_t i = 0; i < GfileFloatKerf.size(); i++) {
    code_lines.push_back(TransformGCodeLine(GfileFloatKerf[i]));
  }
}

std::string Kerf::TransformGCodeLine(const GCodeARRAY_STRUCT &gcode_array) {
  char GCODE[100];
  switch (gcode_array.Name) {
   case G00:
    sprintf(GCODE, "G00 X%.3lf Y%.3lf",
        is_absolute ? gcode_array.X : gcode_array.X - gcode_array.X0,
        is_absolute ? gcode_array.Y : gcode_array.Y - gcode_array.Y0);

    break;
   case G01:
    sprintf(GCODE, "G01 X%.3lf Y%.3lf F%.3lf",
        is_absolute ? gcode_array.X : gcode_array.X - gcode_array.X0,
        is_absolute ? gcode_array.Y : gcode_array.Y - gcode_array.Y0,
        gcode_array.F);

    break;
   case G02:
    sprintf(GCODE, "G02 X%.3lf Y%.3lf I%.3lf J%.3lf F%.3lf",
        is_absolute ? gcode_array.X : gcode_array.X - gcode_array.X0,
        is_absolute ? gcode_array.Y : gcode_array.Y - gcode_array.Y0,
        is_absolute ? gcode_array.I : gcode_array.I - gcode_array.X0,
        is_absolute ? gcode_array.J : gcode_array.J - gcode_array.Y0,
        gcode_array.F);

    break;
   case G03:
    sprintf(GCODE, "G03 X%.3lf Y%.3lf I%.3lf J%.3lf F%.3lf",
        is_absolute ? gcode_array.X : gcode_array.X - gcode_array.X0,
        is_absolute ? gcode_array.Y : gcode_array.Y - gcode_array.Y0,
        is_absolute ? gcode_array.I : gcode_array.I - gcode_array.X0,
        is_absolute ? gcode_array.J : gcode_array.J - gcode_array.Y0,
        gcode_array.F);

    break;
   case G04:
    sprintf(GCODE, "G04 F%.3lf", gcode_array.Delay);
    break;
   case G40:
    sprintf(GCODE, "G40");
    break;
   case G41:
    sprintf(GCODE, "G41");
    break;
   case G42:
    sprintf(GCODE, "G42");
    break;
   case G20:
    sprintf(GCODE, "G20");
    break;
   case G21:
    sprintf(GCODE, "G21");
    break;
   case G90:
    sprintf(GCODE, "G90");
    break;
   case G91:
    sprintf(GCODE, "G91");
    break;
   case G99:
    sprintf(GCODE, "G99 X%.3lf Y%.3lf I%.3lf J%.3lf",
        gcode_array.ScaleFactor, gcode_array.RotateAngle,
        gcode_array.HorizonMirror, gcode_array.VerticalMirror);

    break;
   case M00:
    sprintf(GCODE, "M00");
    break;
   case M02:
    sprintf(GCODE, "M02");
    break;
   case M07:
    sprintf(GCODE, "M07");
    break;
   case M08:
    sprintf(GCODE, "M08");
    break;
   case M09:
    sprintf(GCODE, "M09");
    break;
   case M10:
    sprintf(GCODE, "M10");
    break;
   case M11:
    sprintf(GCODE, "M11");
    break;
   case M12:
    sprintf(GCODE, "M12");
    break;
   default:
    break;
  }
  return std::string(GCODE);
}
