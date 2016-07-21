using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text.RegularExpressions;

namespace SplineyMcSplineFace
{
    public static class Math
    {
        /// <summary>
        /// Factorial function
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static int Fact(this int n)
        {
            int a = 1;
            for (var i = n; i > 0; i--)
            {
                a*=i;
            }
            return a;
        }

        /// <summary>
        /// Factorial function which returns double.
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double FactD(this int n) => (double) Fact(n);


        public static double Pow(this double v, double i) => System.Math.Pow(v, i);

        public static double Bertstein(int i, int n, double u)
        {
            return Bertstein(i,n)(u);
        }

        public static Func<double,double> Bertstein(int i, int n)
        {
            var d = n.FactD()/(i.FactD()*(n - i).FactD());
            return u =>
            {
                Debug.Assert(u >= 0);
                Debug.Assert(u <= 1);
                return d*u.Pow(i)*(1 - u).Pow(n - i);
            };
        }

        public static Func<double,Vector3> Bezier(IReadOnlyList<Vector3> controlPoints)
        {
            var order = controlPoints.Count;
            var degree = order - 1;
            var bersteins = Enumerable.Range(0, order)
                .Select(i => Bertstein(i, degree))
                .ToList();

            return u =>
            {
                Debug.Assert(u >= 0);
                Debug.Assert(u <= 1);
                var sum = Vector3.Zero;
                for (int i = 0; i < order; i++)
                {
                    var ti = (float) bersteins[i](u)*controlPoints[i];
                    sum = sum + ti;
                }
                return sum;

            };
        }

    }

    public class BezierSpline
    {
        public ImmutableList<Vector4> ControlPoints { get; }
        public int Order => ControlPoints.Count;
        public int Degree => Order - 1;

        public BezierSpline(ImmutableList<Vector4> controlPoints)
        {
            ControlPoints = controlPoints;
        }
        public BezierSpline(IEnumerable<Vector3> controlPoints)
        {
            ControlPoints = controlPoints.Select(v=>new Vector4(v,1)).ToImmutableList();
        }

        public BezierSpline Derivative()
        {
            var ctrlPts = new List<Vector4>();
            var degree = Degree;

            for (int i = 0; i < ControlPoints.Count - 1; i++)
            {
                var p = ControlPoints[i+1]-ControlPoints[i];
                ctrlPts.Add(degree * p);
            }
            return new BezierSpline(ctrlPts.ToImmutableList());
        }


        public Vector4 EvaulateAt(double u)
        {
            // deCasteljau algorithm
            var q = ControlPoints.ToList();
            var n = q.Count - 1;
            for (int k = 0; k <= n  ; k++)
            {
                for (int i = 0; i < (n-k) ; i++)
                {
                    q[i] = (float) (1 - u)*q[i] + (float) u*q[i + 1];
                }
                    
            }
            return q[0];
        }

    }

    public class BSpline
    {
        /// <summary>
        /// Return the index of the knot span
        /// </summary>
        /// <param name="k"></param>
        /// <param name="u"></param>
        /// <returns></returns>
        public static int KnotSpanFor(double[] k, double u)
        {
            var lower = -1;
            var i = 0;
            var upper = k.Length;
            do
            {
                i = (upper + lower)/2;
                if (k[i] > u)
                    upper = i;
                else
                {
                    lower = i;
                }
            } while ((upper - lower)>1);

            return lower;
        }


        private class FastDoubleComparer : Comparer<double>
        {
          public override int Compare(double x, double y)
          {
            return x.CompareTo(y);
          }
        }

        /// <summary>
        /// Equation 2.5 in Nurbs book 
        /// </summary>
        /// <param name="k">knot vector</param>
        /// <param name="p">The degree</param>
        /// <param name="i">The ith basic function index</param>
        /// <param name="u">The curve parameter</param>
        /// <returns></returns>
        public static double Basis(double[] k, int p, int i, double u)
        {
            Debug.Assert(i<k.Length);
            Debug.Assert(i>=0);
            if (p == 0)
                return k[i] <= u && u < k[i + 1] ? 1 : 0;


            var a =  (u - k[i])
                   / (k[i + p] - k[i])
                   * Basis(k, p - 1, i, u);

            if (Double.IsNaN(a))
                a = 0;

            var b =  (k[i+p+1] - u)
                   / (k[i + p + 1] - k[i+1])
                   * Basis(k, p - 1, i+1, u);

            if (Double.IsNaN(b))
                b = 0;

            return a + b;
        }

        public static Func<double, double> Basis(double[] k, int p, int i)
        {
            return u => Basis(k, p, i, u);
        }
        public static Func<int, double, double> Basis(double[] k, int p)
        {
            return (i,u) => Basis(k, p, i, u);
        }


    }
}