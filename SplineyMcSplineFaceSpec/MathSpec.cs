using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Linq;
using System.Numerics;
using FluentAssertions;
using MathNet.Numerics;
using SplineyMcSplineFace;
using Xunit;
using Math = SplineyMcSplineFace.Math;

namespace SplineyMcSplineFaceSpec
{
    public class MathSpec
    {
        [Fact]
        public void FactorialShouldWork()
        {
            3.Fact().Should().Be(3*2*1);
            1.Fact().Should().Be(1);
        }

        [Fact]
        public void BernsteinShouldWork()
        {
            Func<double, double> b01 = u => Math.Bertstein(0, 1, u);
            b01(0).Should().Be(1);
            b01(1).Should().Be(0);
        }

        [Fact]
        public void BezierShouldWork()
        {
            var p0 = new Vector3(0,0,0);
            var p1 = new Vector3(1,1,0);
            var p2 = new Vector3(2,0,0);

            var points = new List<Vector3> {p0, p1, p2};

            Math.Bezier(points)(0).Should().Be(p0);
            Math.Bezier(points)(1).Should().Be(p2);
            Math.Bezier(points)(0.5).X.Should().BeGreaterThan(0.9f).And.BeLessThan(1.1f);
            Math.Bezier(points)(0.5).Y.Should().BeGreaterThan(0.4f).And.BeLessThan(0.6f);
        }

        [Fact]
        public void BezierDerivativeShouldWork()
        {
            var p0 = new Vector4(0,0,0,1);
            var p1 = new Vector4(1,1,0,1);
            var p2 = new Vector4(2,0,0,1);

            var list = new[] {p0, p1, p2}.ToImmutableList();

            var spline = new BezierSpline(list);

            var fn = (Func<double,Vector4>) spline.EvaulateAt;

            fn(0).Should().Be(p0);
            fn(1).Should().Be(p2);
            fn(0.5).X.Should().BeGreaterThan(0.9f).And.BeLessThan(1.1f);
            fn(0.5).Y.Should().BeGreaterThan(0.4f).And.BeLessThan(0.6f);



            var dspline = spline.Derivative();

            dspline.Order.Should().Be(spline.Order - 1);

            var dn = (Func<double,Vector4>) dspline.EvaulateAt;
            dn(0).X.Should().BeGreaterThan(0);
            dn(0).Y.Should().BeGreaterThan(0);
            dn(1).X.Should().BeGreaterThan(0);
            dn(1).Y.Should().BeLessThan(0);

        }

        [Fact]
        public void BSplineBasisShouldWork()
        {
            var k = new double[] {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
            BSpline.Basis(k, 0, 0, 0).Should().Be(1);
            BSpline.Basis(k, 0, 0, 0.09).Should().Be(1);
            BSpline.Basis(k, 0, 0, 1).Should().Be(0);

            BSpline.Basis(k, p: 0, i: 4, u: 0.4).Should().Be(1);
            BSpline.Basis(k, p: 0, i: 4, u: 0.49).Should().Be(1);
            BSpline.Basis(k, p: 0, i: 4, u: 0.5).Should().Be(0);

            BSpline.Basis(k, p: 1, i: 0, u: 0).Should().Be(0);
            BSpline.Basis(k, p: 1, i: 0, u: 0.05).Should().BeApproximately(0.5,1e-5);
            BSpline.Basis(k, p: 1, i: 0, u: 0.1).Should().Be(1);
            BSpline.Basis(k, p: 1, i: 0, u: 0.15).Should().BeApproximately(0.5,1e-5);
            BSpline.Basis(k, p: 1, i: 0, u: 0.2).Should().Be(0);

            BSpline.Basis(k, p: 1, i: 4, u: 0.4).Should().Be(0);
            BSpline.Basis(k, p: 1, i: 4, u: 0.5).Should().Be(1);
            BSpline.Basis(k, p: 1, i: 4, u: 0.6).Should().Be(0);


        }


        /// <summary>
        /// Example 2.1 from nurbs book
        /// </summary>
        [Fact]
        public void Example2p1()
        {
            var k = new double[] {0, 0, 0, 1, 1, 1};
            BSpline.Basis(k, i: 0, p: 0).ShouldBe(u=>0,-10,10); // -inf to inf
            BSpline.Basis(k, i: 1, p: 0).ShouldBe(u=>0,-10,-10); 
            BSpline.Basis(k, i: 2, p: 0).ShouldBe(u=>1,0,1); 
            BSpline.Basis(k, i: 3, p: 0).ShouldBe(u=>0,-10,-10); 
            BSpline.Basis(k, i: 4, p: 0).ShouldBe(u=>0,-10,-10); 

            BSpline.Basis(k,i:0, p:1).ShouldBe(u=>0, -10,10);

            BSpline.Basis(k,i:1, p:1).ShouldBe(u=>1-u, 0,1);
            BSpline.Basis(k,i:1, p:1).ShouldBe(u=>0, -10,0);
            BSpline.Basis(k,i:1, p:1).ShouldBe(u=>0, 1,10);

            BSpline.Basis(k,i:0, p:2).ShouldBe(u=>(1-u).Pow(2), 0, 1);
            BSpline.Basis(k,i:1, p:2).ShouldBe(u=>2*u*(1-u), 0, 1);
            BSpline.Basis(k,i:2, p:2).ShouldBe(u=>u*u, 0, 1);

        }

        /// <summary>
        /// Example 2.2 from nurbs book. A selection of the tests
        /// </summary>
        [Fact]
        public void Example2p2()
        {
            var k = new double[] {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
            Func<int, int,Func<double,double>> N = (i, p) => BSpline.Basis(k, i: i, p: p);

            N(0, 0).ShouldBe(u => 0, -10, 10);
            N(1, 0).ShouldBe(u => 0, -10, 10);

            N(2, 0).ShouldBe(u => 1, 0, 1);
            N(3, 0).ShouldBe(u => 1, 1, 2);
            N(4, 0).ShouldBe(u => 1, 2, 3);
            N(5, 0).ShouldBe(u => 1, 3, 4);
            N(6, 0).ShouldBe(u => 0, -10, 10);
            N(7, 0).ShouldBe(u => 1, 4, 5);
            N(8, 0).ShouldBe(u => 0, -10, 10);
            N(9, 0).ShouldBe(u => 0, -10, 10);

            N(0, 1).ShouldBe(u => 0, -10, 10);
            N(1, 1).ShouldBe(u => 1 - u, 0, 1);


            N(2, 1).ShouldBe(u => 0, -10, 0);
            N(2, 1).ShouldBe(u => u, 0, 1);
            N(2, 1).ShouldBe(u => 2 - u, 1, 2);
            N(2, 1).ShouldBe(u => 0, 2, 10);

            N(3, 1).ShouldBe(u => 0, -10, 1);
            N(3, 1).ShouldBe(u => (u - 1), 1, 2);
            N(3, 1).ShouldBe(u => 3 - u, 2, 3);
            N(3, 1).ShouldBe(u => 0, 3, 10);

            N(0,2).ShouldBe(u=>(1-u).Pow(2), 0,1);

            N(1,2).ShouldBe(u=>2*u-3.0/2*u*u, 0,1);
            N(1,2).ShouldBe(u=>0.5*(2-u).Pow(2), 1,2);


            N(5,2).ShouldBe(u=>(u-3).Pow(2),3,4);
            N(5,2).ShouldBe(u=>(5-u).Pow(2),4,5);

        }

        [Fact]
        public void BSplineKnotSearchShortWork()
        {
            {
                var k = new double[] {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
                BSpline.KnotSpanFor(k, 0).Should().Be(2);
                BSpline.KnotSpanFor(k, 0.5).Should().Be(2);
                BSpline.KnotSpanFor(k, 1).Should().Be(3);
                BSpline.KnotSpanFor(k, 2).Should().Be(4);
                BSpline.KnotSpanFor(k, 5).Should().Be(10);
            }
            {
                var k = new double[]{ 0, 0, 1, 1};
                BSpline.KnotSpanFor(k, 0).Should().Be(1);
                BSpline.KnotSpanFor(k, 0.5).Should().Be(1);
                BSpline.KnotSpanFor(k, 1).Should().Be(3);

            }
        }

        [Fact]
        public void PartitionOfUnity()
        {
            var k = new double[] {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
            Func<int, int,Func<double,double>> N = (i, p) => BSpline.Basis(k, i: i, p: p);
            var us = Generate.LinearSpaced(100, 0, 4.9999);
            var p = 2;
            var order = p + 1;
            foreach (var u in us)
            {
                var i = BSpline.KnotSpanFor(k,u);
                Enumerable.Range(i-p, order).Sum(j => N(j, p)(u)).Should().BeApproximately(1.0,1e-5);

            }

        }
    }


    public static class Ext
    {
        public static void ShouldBeAproximatelyEqual(this Vector3 a, Vector3 b, float tol)
        {
            a.X.Should().BeApproximately(b.X, tol);
            a.Y.Should().BeApproximately(b.Y, tol);
            a.Z.Should().BeApproximately(b.Z, tol);
        }

        public enum Bounds

        {
            /// <summary>
            /// Inclusive to inclusive
            /// </summary>
            II,
            /// <summary>
            /// Exclusive to inclusive
            /// </summary>
            EI,
            /// <summary>
            /// Inclusive to exclusive
            /// </summary>
            IE,
            /// <summary>
            /// Exclusive to exclusive
            /// </summary>
            EE

        };

        /// <summary>
        /// A test for comparing two functions over a range of input
        /// </summary>
        /// <param name="fn0"></param>
        /// <param name="fn1"></param>
        /// <param name="lower">lower bound</param>
        /// <param name="upper"></param>
        /// <param name="bounds">Inclusive or Exclusive bounds flags</param>
        /// <param name="eps">Error bound for comparing two function outputs</param>
        public static void ShouldBe
            (this Func<double, double> fn0, Func<double, double> fn1, double lower, double upper, Bounds bounds = Bounds.IE, double eps = 1e-5)
        {
            var steps = 100;
            if (bounds == Bounds.EE || bounds == Bounds.EI)
                lower = lower + (upper - lower)/steps;

            if (bounds == Bounds.EE || bounds == Bounds.IE)
                upper = upper - (upper - lower)/steps;

            var array = Generate.LinearSpaced(steps, lower, upper);
            foreach (var u in array)
            {
                fn0(u).Should().BeApproximately(fn1(u), eps, $"when u = {u}");
            }

        }

        
    }
}