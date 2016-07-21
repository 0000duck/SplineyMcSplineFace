using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Numerics;
using FluentAssertions;
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
        public void BSplineBasisAreSumOfUnity()
        {
            var k = new double[] {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
            BSpline.Basis(k, 0, 0, 0).Should().Be(1);
            BSpline.Basis(k, 0, 0, 0.09).Should().Be(1);
            BSpline.Basis(k, 0, 0, 1).Should().Be(0);

            BSpline.Basis(k, i: 4, p: 0, u: 0.4).Should().Be(1);
            BSpline.Basis(k, i: 4, p: 0, u: 0.49).Should().Be(1);
            BSpline.Basis(k, i: 4, p: 0, u: 0.5).Should().Be(0);

            BSpline.Basis(k, i: 0, p: 1, u: 0).Should().Be(0);
            BSpline.Basis(k, i: 0, p: 1, u: 0.05).Should().BeApproximately(0.5,1e-5);
            BSpline.Basis(k, i: 0, p: 1, u: 0.1).Should().Be(1);
            BSpline.Basis(k, i: 0, p: 1, u: 0.15).Should().BeApproximately(0.5,1e-5);
            BSpline.Basis(k, i: 0, p: 1, u: 0.2).Should().Be(0);

            BSpline.Basis(k, i: 4, p: 1, u: 0.4).Should().Be(0);
            BSpline.Basis(k, i: 4, p: 1, u: 0.5).Should().Be(1);
            BSpline.Basis(k, i: 4, p: 1, u: 0.6).Should().Be(0);


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
        
    }
}