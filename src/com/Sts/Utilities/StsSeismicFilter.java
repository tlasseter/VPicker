package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Filter.*;

import java.nio.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsSeismicFilter
{
	public static boolean UseFocusAGC = true;
	public StsSeismicFilter()
	{
	}

	/** remove DC bias. So that sum doesn't get too large, compute bias for each trace and then average biases for all traces.
	 *  apply same bias to all samples.
	 */
	public static void deBias(float[] fdata, int nCols, int nSlices)
	{
		double sum;
		int start, i;
		sum = 0.0;
		for(int trace = 0, n = 0; trace < nCols; trace++, n += nSlices)
		{
			double traceSum = 0.0;
			for(i = n; i < n + nSlices; i++)
				traceSum += fdata[i];
			traceSum /= nSlices;
			sum += traceSum;
		}
		sum /= nCols;
		for(int n = 0; n < fdata.length; n++)
			fdata[n] -= (float)sum;
	}

	/* sneaky fast box filter */
	public static float[] applyBoxFilter(float[] fdata, int nCols, int nSlices, int filterWindowWidth)
	{
		float[] dest = new float[fdata.length];

		int iptr, itrace, isamp;
		double sum;
		int w2 = filterWindowWidth / 2;
		iptr = 0;
		if (nSlices <= filterWindowWidth) return fdata;
		for(itrace = 0; itrace < nCols; itrace++)
		{
			iptr = itrace * nSlices;
			sum = 0;

			for(int kernel = -w2; kernel <= w2; kernel++)
			{
				sum += fdata[iptr + w2 + kernel];
			}
			dest[iptr + w2] = (float)(sum / filterWindowWidth);
			iptr += w2 + 1;

			for(isamp = w2 + 1; isamp < nSlices - w2; isamp++)
			{
				sum -= fdata[iptr - w2 - 1]; // subtract left sample
				sum += fdata[iptr + w2]; // add right sample

				dest[iptr] = (float)(sum / filterWindowWidth);

				iptr++;
			}
		}
		return dest;
	}

	public static void normalizeAmplitude127(float[] fdata, int nFirstTrace, int nLastTrace, int nSlices)
	{
		float maxAmplitude = 0.0f;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces * nSlices;

		for(int n = 0; n < nFloats; n++)
		{
			float value = fdata[n];
			if(value > maxAmplitude)
				maxAmplitude = value;
			else if(value < -maxAmplitude)
				maxAmplitude = -value;

		}
		for(int n = 0; n < nFloats; n++)
		{
			float value = fdata[n];
			fdata[n] = 127 * (value / maxAmplitude);

		}
	}

	/* This version with bulk get/put, 2.5 times faster in my tests than normalizeAmplitude2, below */

	public static void normalizeAmplitude(FloatBuffer fdata, float[] data, int nFirstTrace, int nLastTrace, int nSlices)
	{
		float maxAmplitude = 0.0f;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces * nSlices;
//              float rawValue = 0, scaledValue = 0;

		fdata.rewind();
		//fdata.get(data);
		for(int n = 0; n < nFloats; n++)
		{
			float value = data[n];
			if(value > maxAmplitude)
				maxAmplitude = value;
			else if(value < -maxAmplitude)
				maxAmplitude = -value;
//                      if(n == 1000) rawValue = value;
		}
		for(int n = 0; n < nFloats; n++)
		{
			float value = data[n];
			data[n] = (value / maxAmplitude);
//                      if(n == 1000) scaledValue = value/maxAmplitude;
		}
		//float checkValue = data.get(1000);
		fdata.rewind();
		fdata.put(data);
		fdata.rewind();
		//data = null;
//              System.out.println("raw " + rawValue + " scaled " + scaledValue + " maxAmp " + maxAmplitude + " check " + checkValue);
	}

	/* in-place max amp normalize */
	public static void normalizeAmplitude(float[] data, int nFirstTrace, int nLastTrace, int nSlices)
	{
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces * nSlices;
        normalizeAmplitude(data, nFloats);
    }

    public static void normalizeAmplitude(float[] data)
    {
        normalizeAmplitude(data, data.length);
    }

    public static void normalizeAmplitude(float[] data, int nFloats)
    {
		float maxAmplitude = 0.0f;
        for(int n = 0; n < nFloats; n++)
		{
			float value = data[n];
			if(value > maxAmplitude)
				maxAmplitude = value;
			else if(value < -maxAmplitude)
				maxAmplitude = -value;
		}
		if(maxAmplitude == 0.0f) return;
//		System.out.println("normalizeAmplitude max " + maxAmplitude + " found at " + maxIndex);
		for(int n = 0; n < nFloats; n++)
			data[n] /= maxAmplitude;
	}

	/* in-place rms  amp normalize */
	public static void normalizeRMSAmplitude(float[] data, int nFirstTrace, int nLastTrace, int nSlices)
	{
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces * nSlices;
		double rms = 0.0;
		float sca = 0;

		for(int n = 0; n < nFloats; n++)
		{
			rms += (data[n] * data[n]);
		}

		rms /= nFloats;
		sca = (float)(1. / Math.sqrt(rms));

		for(int n = 0; n < nFloats; n++)
		{
			data[n] *= sca;
		}
	}

	public static void normalizeAmplitude2(FloatBuffer data, int nFirstTrace, int nLastTrace, int nSlices)
	{
		float maxAmplitude = 0.0f;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces * nSlices;
//		float rawValue = 0, scaledValue = 0;
		data.rewind();
		for(int n = 0; n < nFloats; n++)
		{
			float value = data.get();
			if(value > maxAmplitude)
				maxAmplitude = value;
			else if(value < -maxAmplitude)
				maxAmplitude = -value;
//			if(n == 1000) rawValue = value;
		}
		data.rewind();
		for(int n = 0; n < nFloats; n++)
		{
			float value = data.get(n);
			data.put(value / maxAmplitude);
//			if(n == 1000) scaledValue = value/maxAmplitude;
		}
		float checkValue = data.get(1000);
		data.rewind();
//		System.out.println("raw " + rawValue + " scaled " + scaledValue + " maxAmp " + maxAmplitude + " check " + checkValue);
	}
/*
	public static FloatBuffer applyAGC(FloatBuffer floatBuffer, int nFirstTrace, int nLastTrace, int nCroppedSlices, int windowWidth)
	{
		int i;
		double sum, rms, val;

		if(floatBuffer == null)return null;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		//FloatBuffer agcFloatBuffer = floatBuffer.allocate(nCroppedSlices*nGatherTraces);
		int n = 0;
		float scale = 1.0f;
		float[] data = new float[nGatherTraces * nCroppedSlices];
		float[] agcData = new float[nCroppedSlices];
		floatBuffer.rewind();
		floatBuffer.get(data);
		for(int nTrace = 0; nTrace < nGatherTraces; nTrace++, n += nCroppedSlices)
		{

			// compute initial window for first datum
			sum = 0.0;
			for(i = 0; i < windowWidth + 1; i++)
			{
				//System.out.println("n "+n+ "i "+i+" nCroppedSlices "+nCroppedSlices);
				val = data[n + i];
				sum += val * val;
			}
			int nwin = windowWidth + 1;
			rms = sum / nwin;
			agcData[0] = (rms <= 0.0) ? 0.0f : (float)(scale * data[0] / Math.sqrt(rms));

			// ramping on
			for(i = 1; i <= windowWidth; i++)
			{
				val = data[n + i + windowWidth];
				sum += val * val;
				++nwin;
				rms = sum / nwin;
				agcData[i] = (rms <= 0.0) ? 0.0f : (float)(scale * data[n + i] / Math.sqrt(rms));
			}

			// middle range -- full rms window
			for(i = windowWidth + 1; i <= nCroppedSlices - 1 - windowWidth; i++)
			{
				val = data[n + i + windowWidth];
				sum += val * val;
				val = data[n + i - windowWidth];
				sum -= val * val; // rounding could make sum negative!
				rms = sum / nwin;
				agcData[i] = (rms <= 0.0) ? 0.0f : (float)(scale * data[n + i] / Math.sqrt(rms));
			}

			// ramping off
			for(i = nCroppedSlices - windowWidth; i <= nCroppedSlices - 1; i++)
			{
				val = data[n + i - windowWidth];
				sum -= val * val; // rounding could make sum negative!
				--nwin;
				rms = sum / nwin;
				agcData[i] = (rms <= 0.0) ? 0.0f : (float)(scale * data[n + i] / Math.sqrt(rms));
			}
			System.arraycopy(agcData, 0, data, n, nCroppedSlices);
		}
		floatBuffer.rewind();
		floatBuffer.put(data);
		floatBuffer.rewind();

		return floatBuffer;
	}
*/
	/**
	 * Applies AGC to float[]data.
	 * Uses instantaneousAGC for UseFocusAGC == true, applySuAGC for UseFocusAGC == false.
	 * Converts windowWidthMilliseconds to windowWidth (samples) using ZInc from project.
	 * 
	 */
    public static float[] applyAGC(float[] data, int nFirstTrace, int nLastTrace, int nSlices, int windowWidthMilliseconds, float maxAmp)
    {
        float zInc = StsPreStackLineSetClass.currentProjectPreStackLineSet.zInc;
        int windowWidth = (int) (windowWidthMilliseconds / zInc);
        if (UseFocusAGC)
            return AGC.instantaneousAGC(data, nFirstTrace, nLastTrace, nSlices, windowWidth, maxAmp);
        else
            return applySuAGC(data, nFirstTrace, nLastTrace, nSlices, windowWidth, maxAmp);
    }

    public static float[] applySuAGC(float[] data, int nFirstTrace, int nLastTrace, int nSlices, int windowWidth, float maxAmp)
	{
		int i;
		double sum, rms, val;

		if(data == null)return null;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		//FloatBuffer agcFloatBuffer = floatBuffer.allocate(nCroppedSlices*nGatherTraces);
		int n = 0;
//		float scale = 0.707f;
//		int numFloats = nGatherTraces * nCroppedSlices;
		float minAmpSq = 0.000001f*maxAmp*maxAmp;
	/*
		for (i =0; i < numFloats; i++)
			maxAmp = Math.max(maxAmp,(float)Math.abs(data[i]));
	*/
		// jbw prevent a huge gain from happening based on amp range of data.
//	    maxAmp = .000001f*maxAmp;

		float[] agcData = new float[nSlices];

		for(int nTrace = 0; nTrace < nGatherTraces; nTrace++, n += nSlices)
		{

			/* compute initial window for first datum */
			sum = 0.0;
            for(i = 0; i < windowWidth + 1; i++)
			{
				//System.out.println("n "+n+ "i "+i+" nCroppedSlices "+nCroppedSlices);
				val = data[n + i];
				sum += val * val;
            }
			int nwin = windowWidth + 1;
			rms = sum / nwin;
			agcData[0] = (rms <= minAmpSq) ? data[0] : (float)(data[0] / Math.sqrt(rms));

			/* ramping on */
			for(i = 1; i <= windowWidth; i++)
			{
				val = data[n + i + windowWidth];
				sum += val * val;
				++nwin;
				rms = sum / nwin;
                agcData[i] = (rms <= minAmpSq) ? data[n+i] : (float)(data[n + i] / Math.sqrt(rms));

			}

			/* middle range -- full rms window */
			for(i = windowWidth + 1; i <= nSlices - 1 - windowWidth; i++)
			{
				val = data[n + i + windowWidth];
				sum += val * val;
				val = data[n + i - windowWidth - 1];
				sum -= val * val; /* rounding could make sum negative! */
				rms = sum / nwin;
				agcData[i] = (rms <= minAmpSq) ? data[n+i] : (float)(data[n + i] / Math.sqrt(rms));
			}

			/* ramping off */
			for(i = nSlices - windowWidth; i <= nSlices - 1; i++)
			{
				val = data[n + i - windowWidth - 1];
				sum -= val * val; /* rounding could make sum negative! */
				--nwin;
				rms = sum / nwin;
                agcData[i] = (rms <= minAmpSq) ? data[n+i] : (float)(data[n + i] / Math.sqrt(rms));
			}
			System.arraycopy(agcData, 0, data, n, nSlices);
		}

		return data;
	}
    public static float[] applyAbsAGC(float[] data, int nFirstTrace, int nLastTrace, int nSlices, int windowWidth, float maxAmp)
	{
		int i;
		double sum, rms, val;

		if(data == null)return null;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		//FloatBuffer agcFloatBuffer = floatBuffer.allocate(nCroppedSlices*nGatherTraces);
		int n = 0;
//		float scale = 0.707f;
//		int numFloats = nGatherTraces * nCroppedSlices;
		float minAmp = 0.001f*maxAmp;
	/*
		for (i =0; i < numFloats; i++)
			maxAmp = Math.max(maxAmp,(float)Math.abs(data[i]));
	*/
		// jbw prevent a huge gain from happening based on amp range of data.
//	    maxAmp = .000001f*maxAmp;

		float[] agcData = new float[nSlices];

		for(int nTrace = 0; nTrace < nGatherTraces; nTrace++, n += nSlices)
		{
			/* compute initial window for first datum */
			sum = 0.0;
			for(i = 0; i < windowWidth + 1; i++)
			{
				//System.out.println("n "+n+ "i "+i+" nCroppedSlices "+nCroppedSlices);
				val = Math.abs(data[n + i]);
				sum += val;
			}
			int nwin = windowWidth + 1;
			rms = sum / nwin;
			agcData[0] = (rms <= minAmp) ? data[0] : (float)(data[0] / rms);

			/* ramping on */
			for(i = 1; i <= windowWidth; i++)
			{
				val = Math.abs(data[n + i + windowWidth]);
				sum += val;
				++nwin;
				rms = sum / nwin;
				agcData[i] = (rms <= minAmp) ? data[n+i] : (float)(data[n + i] / rms);
			}

			/* middle range -- full rms window */
			for(i = windowWidth + 1; i <= nSlices - 1 - windowWidth; i++)
			{
				val = Math.abs(data[n + i + windowWidth - 1]);
				sum += val;
				val = Math.abs(data[n + i - windowWidth]);
				sum -= val; /* rounding could make sum negative! */
				rms = sum / nwin;
				agcData[i] = (rms <= minAmp) ? data[n+i] : (float)(data[n + i] / rms);
			}

			/* ramping off */
			for(i = nSlices - windowWidth; i <= nSlices - 1; i++)
			{
				val = Math.abs(data[n + i - windowWidth - 1]);
				sum -= val; /* rounding could make sum negative! */
				--nwin;
				rms = sum / nwin;
				agcData[i] = (rms <= minAmp) ? data[n+i] : (float)(data[n + i] / rms);

			}
		}

		return agcData;
	}

    public static float[] applyAGCtom(float[] inData, int nFirstTrace, int nLastTrace, int nSlices, int windowWidth, float maxAmp)
	{
		int i;
		double sum, rms, val;

		if(inData == null)return null;
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		//FloatBuffer agcFloatBuffer = floatBuffer.allocate(nCroppedSlices*nGatherTraces);
		int n = 0;

	    float[] data = new float[nSlices];
		float[] agcData = new float[nSlices];
		float dataThreshold = 0.001f*maxAmp;
		for(int nTrace = 0; nTrace < nGatherTraces; nTrace++, n += nSlices)
		{
			for(i = 0; i < nSlices; i++)
			{
				float value = inData[n + i];
				if(Math.abs(value) > dataThreshold)
					data[i] = value;
				else
					data[i] = 0.0f;
			}
			// compute initial window for first datum
			sum = 0.0;
			for(i = 0; i < windowWidth + 1; i++)
			{
				//System.out.println("n "+n+ "i "+i+" nCroppedSlices "+nCroppedSlices);
				val = data[i];
				sum += val * val;
			}
			int nwin = windowWidth + 1;
			rms = sum / nwin;
			if(rms > 0.0) agcData[0] = (float)(data[0] / Math.sqrt(rms));
//			agcData[0] = (rms <= 0.001) ? 0.0f : (float)(data[0] / Math.sqrt(rms));

			// ramping on
			for(i = 1; i <= windowWidth; i++)
			{
				val = data[i + windowWidth];
				sum += val * val;
				++nwin;
				rms = sum / nwin;
				if(rms > 0.0) agcData[i] = (float)(data[i] / Math.sqrt(rms));
//				agcData[i] = (rms <= 0.001) ? 0.0f : (float)(data[i] / Math.sqrt(rms));

			}

			// middle range -- full rms window
			for(i = windowWidth + 1; i <= nSlices - 1 - windowWidth; i++)
			{
				val = data[i + windowWidth];
				sum += val * val;
				val = data[i - windowWidth];
				sum -= val * val; // rounding could make sum negative!
				rms = sum / nwin;
				if(rms > 0.0) agcData[i] = (float)(data[i] / Math.sqrt(rms));
//				agcData[i] = (rms <= 0.001) ? 0.0f : (float)(data[i] / Math.sqrt(rms));
			}

			// ramping off
			for(i = nSlices - windowWidth; i <= nSlices - 1; i++)
			{
				val = data[i - windowWidth];
				sum -= val * val; // rounding could make sum negative!
				--nwin;
				rms = sum / nwin;
				if(rms > 0.0) agcData[i] = (float)(data[i] / Math.sqrt(rms));
//				agcData[i] = (rms <= 0.001) ? 0.0f : (float)(data[i] / Math.sqrt(rms));

			}
			System.arraycopy(agcData, 0, inData, n, nSlices);

		}
		return inData;
	}

	public static void butterworthCoefs(float freqlow, float freqhigh, float dt, float[][] b, float[] ynormalizer, int order, boolean[] ift)
	{
		int npole, i, npole2;
		float wdl, wdh, x, wa, s1, sr, si, xn, xd, ang, fnyq;
		float azr, azi, czr, czi, zcr, zci, rtmp, z1r, z1i, z2r, z2i;
		float zt1r, zt1i, zt2r, zt2i, z1mag, z2mag, s1mag, smag;
		final double pi = (double)3.1415926535897932384626433;
		if(order <= 0)return; ;
//  check the validity of the corner frequencies and clean them up if;
//  nessesary.   default to a low pass filter if freqlow is zero and a high;
//  pass filter if freqhigh is greater than the nyquist frequency (fnyq).;
//  otherwise, a pass band filter is generated.;
		ift[0] = false;
		fnyq = 0.5f / dt;
		if(freqlow <= 0.0f)
		{
			freqlow = 0.0f;
			ift[0] = true;
		}
		if(freqhigh > fnyq)freqhigh = fnyq;
		ynormalizer[0] = 1.0f;
		npole = order;
		wdl = (float)pi * freqlow * dt;
		wdh = (float)pi * freqhigh * dt;

		xn = (float)Math.cos(wdl + wdh);
		xd = (float)Math.cos(wdh - wdl);

		x = xn / xd;
		wa = (float)Math.sin(wdh - wdl) / xd;
		s1 = (float)Math.abs(wa);
		npole2 = (npole + 1) / 2;
//  find the filter coefficents.  in the following code, sr and si are related;
//  to the real and imaginary parts of the s plane poles.  there are npole2 pairs;
//  of these poles where the members of each pair are complex conjugates of each;
//  other.  these poles fall on the perimeter of a circle of radius s1.;
		for(i = 1; i <= npole2; i++) //      do 1000 i=1,npole2,1;
		{
			ang = (float)(0.5f * pi * (1.0 + (2.0 * i - 1.0) / npole));
			sr = -s1 * (float)Math.cos(ang);
			si = -s1 * (float)Math.sin(ang);
			azr = 1.0f + sr;
			azi = si;
			czr = 1.0f - sr;
			czi = -si;
			zcr = (azr * czr - azi * czi);
			zci = -(azi * czr + azr * czi);
			zcr = x * x - zcr;
//  find the square root of zc = real(zc) + i*imag(zc) = zcr + i * zci;
			ang = (float)Math.atan2(zci, zcr) / 2.0f;
			rtmp = (float)Math.sqrt(zcr * zcr + zci * zci);
			rtmp = (float)Math.sqrt(rtmp);
			zcr = rtmp * (float)Math.cos(ang);
			zci = rtmp * (float)Math.sin(ang);
			zt1r = x + zcr;
			zt1i = zci;
			zt2r = x - zcr;
			zt2i = -zci;
//  find z1 = z1r + i * z1i = (zt1r + i*zt1i)/(azr + i * azi) and likewise for;
//  z2 = zt2r + i * zt2i over az.;
			rtmp = azr * azr + azi * azi;
			z1r = (zt1r * azr + azi * zt1i) / rtmp;
			z1i = (zt1i * azr - zt1r * azi) / rtmp;
			z2r = (zt2r * azr + azi * zt2i) / rtmp;
			z2i = (zt2i * azr - zt2r * azi) / rtmp;
//  find the filter coeficents for the paired poles located off of the real axis;
			if(i != (npole + 1 - i))
			{
				z1mag = z1r * z1r + z1i * z1i;
				z2mag = z2r * z2r + z2i * z2i;
//            smag          =  sr  * sr  + si  * si;
				smag = s1 * s1;
				s1mag = 1.0f + 2.0f * sr + smag;
				b[0][i - 1] = -2.0f * z1r;
				b[1][i - 1] = z1mag;
				b[0][npole - i] = -2.0f * z2r;
				b[1][npole - i] = z2mag;
				ynormalizer[0] = ynormalizer[0] * smag / s1mag;
			}
			else
			{
//  this section of the routine is aplicable for odd values of order only;
//  and corresponds to the final s plane pole which found on the negative;
//  real axis.;
				rtmp = 1.0f / (1.0f + s1);
				b[0][i - 1] = -2.0f * x * rtmp;
				b[1][i - 1] = (1.0f - s1) * rtmp;
				ynormalizer[0] = ynormalizer[0] * s1 * rtmp;
			}
		}
	}

	public static float[] butterworthFilter(float[] trace, float[][] coefs, float ynorm, int nsect, boolean init, boolean[] ift)
	{
		float xd[] =
			{0.f, 0.f, 0.f};
		int nsamp = trace.length;
		float[][] yd = new float[3][nsect];
		float[] y = new float[trace.length];
		int i = 0;
		int j = 0;
		int jm1 = 0;
		float xtmp = 0.0f;
		float xx0 = 0.0f;
		if((nsect <= 0) || (nsamp <= 0))
		{
			return null;
		}
		if(ift[0])
		{
			xx0 = trace[0];
		}
		else
		{
			xx0 = 0.0f;
		} //  Close else.
// C
// C  Initialize the DELAY array (if appropriate).
// C
		if(init)
		{
			xtmp = trace[0];
			{
				for(i = 0; i < 3; i++)
				{
					xd[i] = xtmp;
					{
						for(j = 0; j < nsect; j++)
						{
							yd[i][j] = 0.0f;
						}
					}
				}
			}
		}
// C
// C  Apply the filter to the trace, X, using a cascade of second order
// C  sections.   Each section takes its input from the previous section
// C  with exception of the specially handled first section.  Note that the
// C  delay array, YD, is updated as soon as its contents are used.
// C
		{
			for(i = 0; i < nsamp; i++)
			{
				xd[0] = trace[(i)];
				yd[0][0] = xd[0] - xd[2] - coefs[0][0] * yd[1][0] - coefs[1][0] * yd[2][0];
				{
					for(j = 1; j < nsect; j++)
					{
						jm1 = j - 1;
						yd[0][j] = yd[0][jm1] - yd[2][jm1]
							- coefs[0][j] * yd[1][j]
							- coefs[1][j] * yd[2][j];

						yd[2][jm1] = yd[1][jm1];
						yd[1][jm1] = yd[0][jm1];
					} //  Close for() loop.
				}
// C
// C Finish updating the delay arrays
// C
				xd[2] = xd[1];
				xd[1] = xd[0];
				yd[2][nsect - 1] = yd[1][nsect - 1];
				yd[1][nsect - 1] = yd[0][nsect - 1];
// C
// C Normalize and store the final, filtered sample
// C
				y[i] = ynorm * yd[0][nsect - 1] + xx0;
			} //  Close for() loop.
		}
// C
		return y;
	}

	public static float[] applyBWFilter(float[] fdata, int nTraces, int nSamps, float bwLow, float bwHigh, int order, float dt)
	{
		float[][] coefs = new float[2][order];

		float[] ynorm = new float[1];
		float[] trace = new float[nSamps];
		boolean[] ift =
			{false};
		butterworthCoefs(bwLow, bwHigh, dt / 1000.f, coefs, ynorm, order, ift);
		for(int itrace = 0; itrace < nTraces; itrace++)
		{
			int n = (nSamps * itrace);
			System.arraycopy(fdata, n, trace, 0, nSamps);
			trace = butterworthFilter(trace, coefs, ynorm[0], order, true, ift);
			System.arraycopy(trace, 0, fdata, n, nSamps);
		}
		return fdata;
	}

/*
    public static float[] applyConvolveFilter(float[] fdata, int nTraces, int nSamps, byte kernel)
    {
        fdata = StsConvolve.convolveFloat2D(fdata, nSamps, nTraces, kernel);
        return fdata;
    }

    public static float[] applyRankFilter(float[] fdata, int nTraces, int nSamps, byte subType, double radius)
    {
        fdata = StsRankFilters.rankFloat(fdata, nSamps, nTraces, radius, subType);
        return fdata;
    }
*/
    public static double [] ricker(double freq, double dt)
	{
	   	double nw = 6./freq/dt;
		int inw =  2 * (int)Math.floor(nw/2.)+1;
		double nc = Math.floor(nw/2.);
		if (nw <= 0) return null;
		double [] w = new double[inw];
		for (int i = 0; i < inw; i++)
		{
			double alpha = (nc -i + 1) * freq * dt * Math.PI;
			double beta = alpha * alpha;
			w[i] =  ((1.-(beta*beta)) * Math.exp(-beta));

		}
		return w;
    }

	public static float[] filter(byte dataType, float[] fdata, int nTraces, int nSamples, float dt, StsFilterProperties filterProperties,
								 StsAGCProperties agcProperties, float dataMin, float dataMax)
	{
		if(fdata == null) return fdata;

        // Box Filter
		if(filterProperties.getApplyBoxFilter(dataType))
			fdata = applyBoxFilter(fdata, nTraces, nSamples, filterProperties.getFilterWindowWidth());

        // Butterworth Filter
        if(filterProperties.getApplyBWFilter(dataType))
			fdata = applyBWFilter(fdata, nTraces, nSamples, filterProperties.getBwLow(),
								  filterProperties.getBWHigh(), filterProperties.getOrder(), dt);

        // Convolution Filter
//        if(filterProperties.getApplyConvolveFilter(dataType))
//            fdata = applyConvolveFilter(fdata, nTraces, nSamples, filterProperties.getKernel());

        // Rank Filter
//        if(filterProperties.getApplyRankFilter(dataType))
 //           fdata = applyRankFilter(fdata, nTraces, nSamples, filterProperties.getSubType(), filterProperties.getFilterRadius());

        // AGC
		if(agcProperties.getApplyAGC(dataType))
		{
			if(dataMin > 0.0 && dataMax > 0.0 || dataMin < 0.0 && dataMax < 0.0) return fdata;
			float maxAmp = Math.min(-dataMin, dataMax);
			fdata = applyAGC(fdata, 0, nTraces-1, nSamples, agcProperties.getWindowWidth(dataType), maxAmp);
//			fdata = applyAGCtom(fdata, 0, nTraces-1, nSamples, agcProperties.getWindowWidth(), maxAmp);
			//normalizeAmplitude(fdata, 0, nTraces-1, nSamples);
		}
//		normalizeAmplitude(fdata, 0, nTraces-1, nSamples);
		return fdata;
	}
/*
    static public void debugCheckAmplitudes(float[] fdata)
	{
		int nValues = fdata.length;
		float min = fdata[0];
		float max = fdata[0];
		for(int n = 1; n < nValues; n++)
		{
			if(fdata[n] < min) min = fdata[n];
			else if(fdata[n] > max) max = fdata[n];
		}
		System.out.println("StsSeismicFilter.debugCheckAmplitudes() min: " + min + " max: " + max);
	}
*/
	public static void main(String[] args)
	{
		float[][] b = new float[2][2];
		float[] ynorm = new float[1];
		boolean[] ift =
			{false};
		butterworthCoefs(20.f, 50.f, 0.002f, b, ynorm, 2, ift);
		System.out.println("ynorm " + ynorm[0]);
		System.out.println(b[0][0] + " " + b[1][0] + " " + b[0][1] + " " + b[1][1]);

		float[] trace = new float[1025];
		float[] result = butterworthFilter(trace, b, ynorm[0], 2, true, ift);
		System.out.println("result " + result);
	}
}
