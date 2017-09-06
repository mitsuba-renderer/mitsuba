/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/testcase.h>

MTS_NAMESPACE_BEGIN

class TestLinearAlgebra : public TestCase {
public:
    MTS_BEGIN_TESTCASE()
    MTS_DECLARE_TEST(test01_basicOperations)
    MTS_DECLARE_TEST(test02_eigenDecomp)
    MTS_DECLARE_TEST(test03_gaussJordan)
    MTS_DECLARE_TEST(test04_lu)
    MTS_DECLARE_TEST(test05_chol)
    MTS_END_TESTCASE()

    void test01_basicOperations() {
        Float refA[] = { 1, 2, 3, 4, 5, 6 };
        Float refB[] = { 7, 8, 9, 10, 11, 12 };
        Float refC[] = { 58, 64, 139, 154 };
        Matrix<2, 3, Float> A(refA);
        Matrix<3, 2, Float> B(refB);
        Matrix<2, 2, Float> C(refC);
        assertEqualsEpsilon(A*B, C, 1e-6);

        Matrix2x2 A2(1, 2, 3, 4), B2(3, 4, 5, 6);
        assertEqualsEpsilon(A2*B2, Matrix2x2(13, 16, 29, 36), 1e-6);
    }

    void test02_eigenDecomp() {
        Matrix4x4 A(
            1.4541, 1.1233, 1.2407, 1.2548,
            1.1233, 0.2597, 0.3819, 1.3917,
            1.2407, 0.3819, 1.1552, 1.1048,
            1.2548, 1.3917, 1.1048, 1.4712
        );
        Matrix4x4 Q;
        Vector4 d;
        A.symEig(Q, (Float *) &d);

        Vector4 refD(-0.823889076095475, 0.130902702868822,
            0.557486242256414, 4.47570013097024);

        Matrix4x4 refQ(
            0.294383137629217, 0.746207030711883, 0.191065628242818, 0.565692108217809,
            -0.789565156329591, 0.139585328270248, -0.459560990857043, 0.381977087957391,
            -0.286639718342894, -0.415237779451033, 0.738328701906588, 0.447533223711729,
            0.455810381691897, -0.501266984714151, -0.45515749947419,  0.577754235510108
        );

        assertEqualsEpsilon(d, refD, 1e-6);
        assertEqualsEpsilon(Q, refQ, 1e-6);
    }

    void test03_gaussJordan() {
        Float orig[5][5] = {
            {0.603220764324061,0.0160774582181293, 0.508210763233408, 0.953606295121746, 0.990510667271985},
            {0.636413405652896, 0.728102544189785, 0.295895331742393, 0.214085933989557, 0.289789452981086},
            {0.757872433764296, 0.176768920514461, 0.222846809207195, 0.437099780859822, 0.285546115075973},
            {0.72645673174838, 0.331824027228289, 0.247810107847795, 0.426238413898285, 0.305370139440192},
            {0.0129043790455411, 0.192371759818767, 0.492432530540112,0.65379652583481, 0.573781874886137}
        };
        Matrix<5, 5, Float> randomMatrix(orig);
        Matrix<5, 5, Float> inverse;
        randomMatrix.invert(inverse);

        Float refInverse[5][5] = {
            {0.0703474910195987,1.40983808579044,5.61251124156711, -5.76223929771622, -0.559883685809967},
            {-0.0613980284171623, -1.63235458842233, -10.5399802908894,12.4851221198693, -0.468956103463839},
            {-2.29568951232449,7.39978908961836,27.3281328010324, -33.1622231774198, 4.27485563386439},
            {-0.885020418512179, -6.82656789249295, -15.71882942499, 23.104992960391, 0.501523848970031},
            {2.99765294474105, 1.9434481868496, -2.13526601478719, -1.92274014105124, -2.32759821076299}
        };
        Matrix<5, 5, Float> refInverseMatrix(refInverse);
        assertEqualsEpsilon(inverse, refInverseMatrix, 1e-4);
    }

    void test04_lu() {
        Float orig[5][5] = {
            {0.603220764324061,0.0160774582181293, 0.508210763233408, 0.953606295121746, 0.990510667271985},
            {0.636413405652896, 0.728102544189785, 0.295895331742393, 0.214085933989557, 0.289789452981086},
            {0.757872433764296, 0.176768920514461, 0.222846809207195, 0.437099780859822, 0.285546115075973},
            {0.72645673174838, 0.331824027228289, 0.247810107847795, 0.426238413898285, 0.305370139440192},
            {0.0129043790455411, 0.192371759818767, 0.492432530540112,0.65379652583481, 0.573781874886137}
        };
        Matrix<5, 5, Float> randomMatrix(orig);
        Matrix<5, 5, Float> lu;
        int piv[5], pivsign;
        if (!randomMatrix.lu(lu, piv, pivsign))
            Log(EError, "Could not generate LU decomposition!");

        Float det = lu.luDet(pivsign);
        assertEqualsEpsilon(det, 0.00294638676819103, 1e-6);
        assertEqualsEpsilon(randomMatrix.det(), 0.00294638676819103, 1e-6);

        Float b[5] = { 0.088842801947558, 0.717330636910085,
            0.689651610365336, 0.76615421808756, 0.274406750300772 };
        Float xref[5] = { 0.319847529827911, 0.991537144805428,
            -0.283297199558031, 2.02355933108543, -1.92399896347868 };

        Matrix<5, 1, Float> B(b);
        Matrix<5, 1, Float> X;
        lu.luSolve(B, X, piv);

        Matrix<5, 1, Float> Xref(xref);
        assertEqualsEpsilon(X, Xref, 1e-5);
    }

    void test05_chol() {
        Float orig[5][5] = {
            {5.51302956379742,     -1.36613599367125,     -0.56818386886837,      1.43247259160815,    -0.341242075046427},
            {-1.36613599367125,      4.32377082646249,     -1.78941813087798,   -0.0389694848978303,     -1.36129441356894},
            {-0.56818386886837,     -1.78941813087798,      3.16179437285643,    -0.169145503847853,    -0.175824617890457},
            {1.43247259160815,   -0.0389694848978303,    -0.169145503847853,      3.43678526823342,   -0.0421265728657415},
            {-0.341242075046427,     -1.36129441356894,    -0.175824617890457,   -0.0421265728657415,      2.75055078833228}
        };
        Float cholRef[5][5] = {
            {2.34798414896639,                     0,                     0,                     0,                     0},
            {-0.581833567433853,      1.99630672149089,                     0,                     0,                     0},
            {-0.241987949159917,    -0.966892923734769,      1.47253328632987,                     0,                     0},
            {0.610086142293059,     0.158291863826668,    0.0893285569680872,      1.74113303971054,                    0},
            {-0.145334062496395,    -0.724264780584955,    -0.618852021434583,     0.124324850430034,      1.34403676785521}
        };
        Matrix<5, 5, Float> M(orig), L;
        if (!M.chol(L))
            Log(EError, "Could not generate Cholesky decomposition!");
        Matrix<5, 5, Float> Lref(cholRef);
        assertEqualsEpsilon(L, Lref, 1e-6);
        assertEqualsEpsilon(L.cholDet(), 260.892331265789, 1e-4);

        Float b[5] = { 0.331003261288103, 0.601514197755584, 0.0474150216556957, 0.363418224633018, 0.47712128258488 };
        Float xref[5] = { 0.25931313190664, 0.544146878658952, 0.398994177053345, 0.0296072370438656, 0.500901199816733 };

        Matrix<5, 1, Float> B(b);
        Matrix<5, 1, Float> X;
        L.cholSolve(B, X);

        Matrix<5, 1, Float> Xref(xref);
        assertEqualsEpsilon(X, Xref, 1e-5);
    }
};

MTS_EXPORT_TESTCASE(TestLinearAlgebra, "Testcase for Linear Algebra routines")
MTS_NAMESPACE_END
