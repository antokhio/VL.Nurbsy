/*
 * Ported by:
 * 2025 - Anton Kalabuhov (antokhio)
 *
 * Original Author:
 * 2024 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 * Based on LNLib: https://github.com/BIMCoderLiang/LNLib
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

namespace Nurbsy.Algorithm
{
    public delegate double IntegrationFunction(double t, object customData);
    public delegate double BinaryIntegrationFunction(double u, double v, object customData);

    public static class Integrator
    {
        public static double Simpson(
            IntegrationFunction function,
            object customData,
            double start,
            double end
        )
        {
            double st = function(start, customData);
            double mt = function((start + end) / 2.0, customData);
            double et = function(end, customData);
            double result = ((end - start) / 6.0) * (st + 4 * mt + et);
            return result;
        }

        public static double Simpson(
            BinaryIntegrationFunction function,
            object customData,
            double uStart,
            double uEnd,
            double vStart,
            double vEnd
        )
        {
            double du = uEnd - uStart;
            double dv = vEnd - vStart;
            double hdu = 0.5 * du;
            double hdv = 0.5 * dv;

            // Sample 9 points with weights
            int sampleNumber = 9;
            int patches = 4;
            // uvw: u, v, w (weight) flattened
            ReadOnlySpan<double> uvw =
                stackalloc double[27] {
                    uStart,
                    vStart,
                    1.0,
                    uStart,
                    vStart + hdv,
                    4.0,
                    uStart,
                    vEnd,
                    1.0,
                    uStart + hdu,
                    vStart,
                    4.0,
                    uStart + hdu,
                    vStart + hdv,
                    16.0,
                    uStart + hdu,
                    vEnd,
                    4.0,
                    uEnd,
                    vStart,
                    1.0,
                    uEnd,
                    vStart + hdv,
                    4.0,
                    uEnd,
                    vEnd,
                    1.0,
                };

            double sum = 0;
            for (int i = 0; i < sampleNumber; ++i)
            {
                int baseIdx = i * 3;
                double u = uvw[baseIdx];
                double v = uvw[baseIdx + 1];
                double w = uvw[baseIdx + 2];
                double f = function(u, v, customData);
                sum += w * f;
            }

            sum *= du * dv / (sampleNumber * patches);
            return sum;
        }

        public static readonly IReadOnlyList<double> GaussLegendreAbscissae = new[]
        {
            -0.0640568928626056260850430826247450385909,
            0.0640568928626056260850430826247450385909,
            -0.1911188674736163091586398207570696318404,
            0.1911188674736163091586398207570696318404,
            -0.3150426796961633743867932913198102407864,
            0.3150426796961633743867932913198102407864,
            -0.4337935076260451384870842319133497124524,
            0.4337935076260451384870842319133497124524,
            -0.5454214713888395356583756172183723700107,
            0.5454214713888395356583756172183723700107,
            -0.6480936519369755692524957869107476266696,
            0.6480936519369755692524957869107476266696,
            -0.7401241915785543642438281030999784255232,
            0.7401241915785543642438281030999784255232,
            -0.8200019859739029219539498726697452080761,
            0.8200019859739029219539498726697452080761,
            -0.8864155270044010342131543419821967550873,
            0.8864155270044010342131543419821967550873,
            -0.9382745520027327585236490017087214496548,
            0.9382745520027327585236490017087214496548,
            -0.9747285559713094981983919930081690617411,
            0.9747285559713094981983919930081690617411,
            -0.9951872199970213601799974097007368118745,
            0.9951872199970213601799974097007368118745,
        };

        public static readonly IReadOnlyList<double> GaussLegendreWeights = new[]
        {
            0.1279381953467521569740561652246953718517,
            0.1279381953467521569740561652246953718517,
            0.1258374563468282961213753825111836887264,
            0.1258374563468282961213753825111836887264,
            0.121670472927803391204463153476262425607,
            0.121670472927803391204463153476262425607,
            0.1155056680537256013533444839067835598622,
            0.1155056680537256013533444839067835598622,
            0.1074442701159656347825773424466062227946,
            0.1074442701159656347825773424466062227946,
            0.0976186521041138882698806644642471544279,
            0.0976186521041138882698806644642471544279,
            0.086190161531953275917185202983742667185,
            0.086190161531953275917185202983742667185,
            0.0733464814110803057340336152531165181193,
            0.0733464814110803057340336152531165181193,
            0.0592985849154367807463677585001085845412,
            0.0592985849154367807463677585001085845412,
            0.0442774388174198061686027482113382288593,
            0.0442774388174198061686027482113382288593,
            0.0285313886289336631813078159518782864491,
            0.0285313886289336631813078159518782864491,
            0.0123412297999871995468056670700372915759,
            0.0123412297999871995468056670700372915759,
        };

        public static double[] ChebyshevSeries(int size)
        {
            var series = new double[size];

            int lenw = size - 1;
            int j,
                k,
                l,
                m;
            double cos2,
                sin1,
                sin2,
                hl;

            cos2 = 0;
            sin1 = 1;
            sin2 = 1;
            hl = 0.5;
            k = lenw;
            l = 2;
            while (l < k - l - 1)
            {
                series[0] = hl * 0.5;
                for (j = 1; j <= l; j++)
                {
                    series[j] = hl / (1 - 4 * j * j);
                }
                series[l] *= 0.5;
                FFT.Dfct(l, 0.5 * cos2, sin1, series);
                cos2 = Math.Sqrt(2 + cos2);
                sin1 /= cos2;
                sin2 /= 2 + cos2;
                series[k] = sin2;
                series[k - 1] = series[0];
                series[k - 2] = series[l];
                k -= 3;
                m = l;
                while (m > 1)
                {
                    m >>= 1;
                    for (j = m; j <= l - m; j += (m << 1))
                    {
                        series[k] = series[j];
                        k--;
                    }
                }
                hl *= 0.5;
                l *= 2;
            }
            return series;
        }

        public static double ClenshawCurtisQuadrature(
            IntegrationFunction function,
            object customData,
            double start,
            double end,
            Span<double> series,
            double epsilon
        )
        {
            double integration;
            int j,
                k,
                l;
            double err,
                esf,
                eref,
                erefh,
                hh,
                ir,
                iback,
                irback,
                ba,
                ss,
                x,
                y,
                fx,
                errir;
            int lenw = series.Length - 1;
            esf = 10;
            ba = 0.5 * (end - start);
            ss = 2 * series[lenw];
            x = ba * series[lenw];
            series[0] = 0.5 * function(start, customData);
            series[3] = 0.5 * function(end, customData);
            series[2] = function(start + x, customData);
            series[4] = function(end - x, customData);
            series[1] = function(start + ba, customData);
            eref =
                0.5
                * (
                    Math.Abs(series[0])
                    + Math.Abs(series[1])
                    + Math.Abs(series[2])
                    + Math.Abs(series[3])
                    + Math.Abs(series[4])
                );
            series[0] += series[3];
            series[2] += series[4];
            ir = series[0] + series[1] + series[2];
            integration =
                series[0] * series[lenw - 1]
                + series[1] * series[lenw - 2]
                + series[2] * series[lenw - 3];
            erefh = eref * Math.Sqrt(epsilon);
            eref *= epsilon;
            hh = 0.25;
            l = 2;
            k = lenw - 5;
            do
            {
                iback = integration;
                irback = ir;
                x = ba * series[k + 1];
                y = 0;
                integration = series[0] * series[k];
                for (j = 1; j <= l; j++)
                {
                    x += y;
                    y += ss * (ba - x);
                    fx = function(start + x, customData) + function(end - x, customData);
                    ir += fx;
                    integration += series[j] * series[k - j] + fx * series[k - j - l];
                    series[j + l] = fx;
                }
                ss = 2 * series[k + 1];
                err = esf * l * Math.Abs(integration - iback);
                hh *= 0.25;
                errir = hh * Math.Abs(ir - 2 * irback);
                l *= 2;
                k -= l + 2;
            } while ((err > erefh || errir > eref) && k > 4 * l);
            integration *= end - start;
            if (err > erefh || errir > eref)
            {
                err *= -Math.Abs(end - start);
            }
            else
            {
                err = eref * Math.Abs(end - start);
            }
            return integration;
        }

        // Just an overload that accepts array/list implicitly as Span or copies logic?
        // The original logic is identical for Quadrature2, just passed vector by value in C++ copy which was weird.
        // Assuming we keep one robust implementation with Span.
    }
}
