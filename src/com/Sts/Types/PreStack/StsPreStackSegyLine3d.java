package com.Sts.Types.PreStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.StsSEGYFormat.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.util.*;

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
public class StsPreStackSegyLine3d extends StsPreStackSegyLine
{
  /** x of first CDP for a straight 3d line */
    public float firstX;
    /** y of first CDP for a straight 3d line */
    public float firstY;
    /** x of last CDP for a straight 3d line  */
    public float lastX;
    /** y of last CDP for a straight 3d line  */
    public float lastY;

    transient int nOutputSamplesPerWrite = 0;

    public String getGroupname()
    {
        return group3dPrestack;
    }

    public StsPreStackSegyLine3d()
	{
		super();
	}

	public StsPreStackSegyLine3d(boolean persistent)
	{
		super(persistent);
	}

	protected StsPreStackSegyLine3d(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws FileNotFoundException, StsException
	{
		super(seismicWizard, file, stsDirectory, frame, segyFormat);
	}

	static public StsPreStackSegyLine3d constructor(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat)
	{
		try
		{
			return new StsPreStackSegyLine3d(seismicWizard, file, stsDirectory, frame, segyFormat);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyVolume.constructTraceAnalyzer(file, stsDirectory) failed.", e,
										 StsException.WARNING);
			return null;
		}
	}

    public String createFilename(String format)
	{
        String lineName;
        if(isInline)
            lineName = "I" + (int)rowNumMin + "X" + (int)colNumMin + "-" + (int)colNumMax;
        else
            lineName = "X" + (int)colNumMin + "I" + (int)rowNumMin + "-" + (int)rowNumMax;

        String fullStemname = stemname + "." + "line." + lineName;
        return createFilename(format, fullStemname);
	}

    /** computes volume geometry from trace CDP ilines, crossline, X, and Y values */
	public boolean analyzeGeometry()
	{
		try
		{
            originSet = false;
        /*
            setBytesPerSample(StsSegyIO.getBytesPerSample(getSampleFormat()));
			setTraceHeaderSize(getSegyFormat().traceHeaderSize);
			setTraceHeader(new byte[traceHeaderSize]);
			setBytesPerTrace(nCroppedSlices * getBytesPerSample() + traceHeaderSize);
			setSegyFileSize(randomAccessSegyFile.length());
	    */
//            float nTotalTracesFloat = (float)((float)(segyFileSize - fileHeaderSize) / (float)bytesPerTrace);
//			setNTotalTraces((int)(((segyFileSize - fileHeaderSize)/numVolumes) / bytesPerTrace));
//            if(nTotalTracesFloat != nTotalTraces)
//            {
//                new StsMessage(frame, StsMessage.WARNING, "Unexpected file size, please verify header sizes.");
//                return false;
//            }
//			nTotalTraces = 10555;
//			if(getNTotalTraces() == 0)
//				return false;
			/*
			 int nHeaderTraces = (int)segyFormat.getBinaryRec("NTRACES").getHdrValue(getBinaryHeader(), getIsLittleEndian);
			   if(nHeaderTraces != nTotalTraces)
			   {
			 new StsMessage(frame, StsMessage.WARNING, "Number of traces found: " + nTotalTraces + " doesn't agree with binary header: " + nHeaderTraces);
			   }
			 */
//            tMin = segyFormat.startT * verticalScalar;

			// Get first trace and first trace in next gather to determine if this is line or xline
			firstTraceHeader = getTrace(0);
			if(firstTraceHeader == null)return false;

//			double dx = secondTraceHeader.x - firstTraceHeader.x;
//			double dy = secondTraceHeader.y - firstTraceHeader.y;
//			xInc = (float)Math.sqrt(dx*dx + dy*dy);

			lastTraceHeader = getTrace(getNTotalTraces() - 1);
			if(lastTraceHeader == null)return false;

			xOrigin = firstTraceHeader.x;
			yOrigin = firstTraceHeader.y;

			rowNumMin = firstTraceHeader.iLine;
			colNumMin = firstTraceHeader.xLine;
			rowNumMax = lastTraceHeader.iLine;
			colNumMax = lastTraceHeader.xLine;

			secondTraceHeader = findFirstTraceSecondGather();
			if(secondTraceHeader == null)return false;

//            if(((rowNumMax - rowNumMin) == 0) || ((colNumMax - colNumMin) == 0)) return false;
			findLastTraceFirstLine();
			findFirstTraceLastLine();

//			isTracesPerGatherEqual(nTracesPerGather, nTotalTraces);

			if(!analyzeAngle(firstTraceHeader, firstLineLastTraceHeader, lastLineFirstTraceHeader))
			{
				System.out.println("Failed to analyze angle....");
				return false;
			}
            setLocalXYRange();
            return true;
		}
		// again just use (Exception e) and you get everything
		catch(Exception e)
		{
			StsMessageFiles.errorMessage("Couldn't analyze geometry: select a different number of bytes per sample.");
			return false;
		}
	}

	private StsSEGYFormat.TraceHeader findFirstTraceSecondGather()
	{
		secondTraceHeader = null;

		rowNumMin = firstTraceHeader.iLine;
		colNumMin = firstTraceHeader.xLine;
		nTracesPerGather = 1;
		while(nTracesPerGather < 1000)
		{
			secondTraceHeader = getTrace(nTracesPerGather);
			if(secondTraceHeader.iLine == rowNumMin && secondTraceHeader.xLine != colNumMin)
			{
				isInline = true;
//				cdpNumMin = colNumMin;
				colNumInc = secondTraceHeader.xLine - colNumMin;
				break;
			}
			if(secondTraceHeader.xLine == colNumMin && secondTraceHeader.iLine != rowNumMin)
			{
				isInline = false;
//				cdpNumMin = rowNumMin;
				rowNumInc = secondTraceHeader.iLine - rowNumMin;
				break;
			}
			secondTraceHeader = null;
			nTracesPerGather++;
		}
		if(secondTraceHeader == null)
		{
			StsMessageFiles.errorMessage("Unable to find next gather.");
			return null;
		}
		return secondTraceHeader;
	}

	protected void findLastTraceFirstLine()
	{
		// binary increment search to find it; begin by assuming the lines numbers are spaced by one
		// so we can estimate the number of lines
		int numberOfLines;
		if(isInline)
			numberOfLines = (int)Math.abs(rowNumMax - rowNumMin);
		else
			numberOfLines = (int)Math.abs(colNumMax - colNumMin);
		numberOfLines = StsMath.minMax(numberOfLines, 1, 1000);
		if(numberOfLines == 1)
		{
			firstLineLastTraceHeader = lastTraceHeader;
			return;
		}
		int estNTracesPerLine = getNTotalTraces() / numberOfLines;
		int step = estNTracesPerLine / 10;
		int nTrace = estNTracesPerLine - step / 2;
		int iter = 0;
		StsSEGYFormat.TraceHeader prevTraceHeader = null;
		StsSEGYFormat.TraceHeader nextTraceHeader = getTrace(nTrace);
		float firstLineNum;
		if(isInline)
		{
			firstLineNum = firstTraceHeader.iLine;
			while(step != 0)
			{
				prevTraceHeader = nextTraceHeader;
				nTrace += step;
				nextTraceHeader = getTrace(nTrace);
				step = searchForEndOfLine(firstLineNum, prevTraceHeader.iLine, nextTraceHeader.iLine, step);
			}
			firstLineLastTraceHeader = prevTraceHeader;
			rowNumInc = nextTraceHeader.iLine - prevTraceHeader.iLine;
//			double dx = nextTraceHeader.x - firstTraceHeader.x;
//			double dy = nextTraceHeader.y - firstTraceHeader.y;
//			yInc = (float)Math.sqrt(dx*dx + dy*dy);
		}
		else
		{
			firstLineNum = (float)firstTraceHeader.xLine;
			while(step != 0)
			{
				prevTraceHeader = nextTraceHeader;
				nTrace += step;
				nextTraceHeader = getTrace(nTrace);
				step = searchForEndOfLine(firstLineNum, prevTraceHeader.xLine, nextTraceHeader.xLine, step);
			}
			firstLineLastTraceHeader = prevTraceHeader;
			colNumInc = nextTraceHeader.xLine - prevTraceHeader.xLine;
		}
	}

	private int searchForEndOfLine(float firstLineNum, float prevLineNum, float nextLineNum, int step)
	{
		// if both prev and next are on same line as firstLine, search forward; if step is negative, flip & halve
		// if both prev and next are on same line, but not on firstLine, search backward; if step is positive, flip & halve
		// if prev on firstLine and next is not, search backwards; if step is postive, flip & halve
		// if prev not on firstLine and next is, search forward; if step is negative, flip & halve

		boolean isPrev = firstLineNum == prevLineNum;
		boolean isNext = firstLineNum == nextLineNum;
		if(isPrev)
		{
			if(isNext)
			{
				if(step < 0)return flipStep(step);
				else return step;
			}
			else
			{
				if(step == 1)return 0;
				if(step > 0)return flipStep(step);
				else return step;
			}
		}
		else // !isPrev
		{
			if(isNext)
			{
				if(step < 0)return flipStep(step);
				else return step;
			}
			else
			{
				if(step > 0)return flipStep(step);
				else return step;
			}
		}
	}

	private int flipStep(int step)
	{
		if(step == 1)return -1;
		else if(step == -1)return 1;
		else return -step / 2;
	}

	protected void findFirstTraceLastLine()
	{
		// binary increment search to find it; begin by assuming the lines numbers are spaced by one
		// so we can estimate the number of lines
		int numberOfLines;
		if(isInline)
			numberOfLines = (int)Math.abs(rowNumMax - rowNumMin);
		else
			numberOfLines = (int)Math.abs(colNumMax - colNumMin);
		numberOfLines = StsMath.minMax(numberOfLines, 1, 1000);
		if(numberOfLines <= 1)
		{
			lastLineFirstTraceHeader = null;
			return;
		}

		int estNTracesPerLine = getNTotalTraces() / numberOfLines;
		int step = estNTracesPerLine / 10;
		int nTrace = getNTotalTraces() - estNTracesPerLine - step / 2;
		StsSEGYFormat.TraceHeader prevTraceHeader = null;
		StsSEGYFormat.TraceHeader nextTraceHeader = getTrace(nTrace);
		float lastLineNum;
		if(isInline)
		{
			lastLineNum = lastTraceHeader.iLine;
			while(step != 0)
			{
				prevTraceHeader = nextTraceHeader;
				int nextNTrace = nTrace + step;
				if(step > 0 && nextNTrace >= getNTotalTraces())
				{
					nextNTrace = (nTrace + getNTotalTraces()) / 2;
					step = nextNTrace - nTrace;
				}
				else if(step < 0 && nextNTrace < 0)
				{
					nextNTrace = nTrace / 2;
					step = nTrace - nextNTrace;
				}
				nTrace = nextNTrace;
				nextTraceHeader = getTrace(nTrace);
				step = searchForStartOfLine(lastLineNum, prevTraceHeader.iLine, nextTraceHeader.iLine, step);
			}
			lastLineFirstTraceHeader = nextTraceHeader;
		}
		else
		{
			lastLineNum = (float)firstTraceHeader.xLine;
			while(step != 0)
			{
				prevTraceHeader = nextTraceHeader;
				int nextNTrace = nTrace + step;
				if(step > 0 && nextNTrace >= getNTotalTraces())
				{
					nextNTrace = (nTrace + getNTotalTraces()) / 2;
					step = nextNTrace - nTrace;
				}
				else if(step < 0 && nextNTrace < 0)
				{
					nextNTrace = nTrace / 2;
					step = nTrace - nextNTrace;
				}
				nTrace = nextNTrace;
				nextTraceHeader = getTrace(nTrace);
				step = searchForStartOfLine(lastLineNum, prevTraceHeader.xLine, nextTraceHeader.xLine, step);
			}
			lastLineFirstTraceHeader = nextTraceHeader;
		}
	}

	private int searchForStartOfLine(float lastLineNum, float prevLineNum, float nextLineNum, int step)
	{
		// if both prev and next are on same line as lastLine, search backward; if step is positive, flip & halve
		// if both prev and next are on same line, but not on firstLine, search forward; if step is negative, flip & halve
		// if prev on lastLine and next is not, search forwards; if step is negative, flip & halve
		// if prev not on lastLine and next is, search backward; if step is positive, flip & halve


		boolean isPrev = lastLineNum == prevLineNum;
		boolean isNext = lastLineNum == nextLineNum;
		if(isPrev)
		{
			if(isNext)
			{
				if(step > 0)return flipStep(step);
				else return step;
			}
			else
			{

				if(step < 0)return flipStep(step);
				else return step;
			}
		}
		else // !isPrev
		{
			if(isNext)
			{
				if(step == 1)return 0;
				if(step > 0)return flipStep(step);
				else return step;
			}
			else
			{
				if(step < 0)return flipStep(step);
				else return step;
			}
		}
	}

	/** Given an origin trace (arbitrary) and two other traces, compute the line and xline spacings
	 *  and rotation angle of grid.  Any 3 traces can by supplied as long as they are not in a straight line.
	 *  The line and xline numbering increments already need to have been computed.
	 *  We assume lines and xlines are orthogonal. We have two vectors from origin to each of the two points.
	 *  Each vector has components in the line and xline directions in which we know the index difference in
	 *  each component but not the index spacing.  We have 3 unknowns: yInc, xInc, and lineAngle
	 *  (the angle from the +X axis to the +line direction.  Using Pythagorean theorem, we solve two equations
	 *  for the two spacing unknowns.  The angle is then computed using the dot and cross products between
	 *  the same vector in the rotated and unrotated coordinate systems.
	 */
	public boolean analyzeAngle(StsSEGYFormat.TraceHeader originTrace, StsSEGYFormat.TraceHeader traceA, StsSEGYFormat.TraceHeader traceB)
	{
		double dLineA = 0.0, dXLineA = 0.0;
		double dLineB = 0.0, dXLineB = 0.0;
		double dLineSqA = 0.0, dXLineSqA = 0.0;
		double dLineSqB = 0.0, dXLineSqB = 0.0;
		double lengthSqA = 0.0, lengthSqB = 0.0;
		double dxA = 0.0, dyA = 0.0, dxB = 0.0, dyB = 0.0;
		boolean traceAok = false, traceBok = false;

		if(traceA != null)
		{
			if(rowNumInc != 0.0f)
			{
				dLineA = (traceA.iLine - originTrace.iLine) / rowNumInc;
				dLineSqA = dLineA * dLineA;
			}
			if(colNumInc != 0.0f)
			{
				dXLineA = (traceA.xLine - originTrace.xLine) / colNumInc;
				dXLineSqA = dXLineA * dXLineA;
			}
			dxA = traceA.x - originTrace.x;
			dyA = traceA.y - originTrace.y;
			lengthSqA = dxA * dxA + dyA * dyA;
			if(lengthSqA > 1.0f)traceAok = true;
		}
		if(traceB != null)
		{
			if(rowNumInc != 0.0f)
			{
				dLineB = (traceB.iLine - originTrace.iLine) / rowNumInc;
				dLineSqB = dLineB * dLineB;
			}
			if(colNumInc != 0.0f)
			{
				dXLineB = (traceB.xLine - originTrace.xLine) / colNumInc;
				dXLineSqB = dXLineB * dXLineB;
			}
			dxB = traceB.x - originTrace.x;
			dyB = traceB.y - originTrace.y;
			lengthSqB = dxB * dxB + dyB * dyB;
			if(lengthSqB > 1.0f)traceBok = true;
		}

		if(!traceAok && !traceBok)return false;
		if(traceAok && traceBok)
		{
			// check if angle between two vectors is too close ( less than 1 degree)
			double crossIndexes = dXLineA * dLineB - dLineA * dXLineB;
			double dotIndexes = dXLineA * dXLineB + dLineA * dLineB;
			double indexAngle = StsMath.atan2(dotIndexes, crossIndexes);
			if(Math.abs(indexAngle) < 1)return false;

			double crossSq = dXLineSqA * dLineSqB - dXLineSqB * dLineSqA;

			//        if(Math.abs(crossSq) < 0.1*lengthSqA) return false;

			double yIncSq = (lengthSqB * dXLineSqA - lengthSqA * dXLineSqB) / crossSq;
			double yInc = Math.sqrt(yIncSq) * horizontalScalar;

			double xIncSq;
			if(dXLineSqA > dXLineSqB)
				xIncSq = (lengthSqA - yIncSq * dLineSqA) / dXLineSqA;
			else
				xIncSq = (lengthSqB - yIncSq * dLineSqB) / dXLineSqB;
			double xInc = Math.sqrt(xIncSq) * horizontalScalar;

			// check if xLines are 90 deg CCW from lines; if not change sign of yInc
			double crossLines = dxA * dyB - dyA * dxB;

			isXLineCCW = (crossIndexes * crossLines >= 0);
//			if(!isXLineCCW) flipRowNumOrder();

			// compute angle from same vector in rotated coordinate system to vector in unrotated coordinate system
			// this is the rotation angle from +X in unrotated coordinates to +Line direction
			double cosL = (xInc * dXLineA * dxA + yInc * dLineA * dyA);
			double sinL = (xInc * dXLineA * dyA - yInc * dLineA * dxA);
			double lineAngle = StsMath.atan2(cosL, sinL);
			this.yInc = (float)yInc;
			this.xInc = (float)xInc;
			this.angle = (float)lineAngle;
		}
		else if(traceAok)
		{
			double lengthA = Math.sqrt(lengthSqA);
			if(dLineA != 0.0)
			{
				this.yInc = (float)(lengthA / dLineA);
				this.xInc = 0.0f;
				this.angle = (float)StsMath.atan2(dxA, dyA);
			}
			else if(dXLineA != 0.0)
			{
				this.yInc = 0.0f;
				this.xInc = (float)(lengthA / dXLineA);
				this.angle = (float)StsMath.atan2(dxA, dyA);
			}
			else
				return false;
		}
		else if(traceBok)
		{
			double lengthB = Math.sqrt(lengthSqB);
			if(dLineB != 0.0)
			{
				this.yInc = (float)(lengthB / dLineB);
				this.xInc = 0.0f;
				this.angle = (float)StsMath.atan2(dxB, dyB);
			}
			else if(dXLineB != 0.0)
			{
				this.yInc = 0.0f;
				this.xInc = (float)(lengthB / dXLineB);
				this.angle = (float)StsMath.atan2(dxB, dyB);
			}
			else
				return false;
		}
		return true;
	}

	// JKF 27JUNE2006 - resurrected this
	protected void setLocalXYRange()
	{
        float firstX = firstTraceHeader.x;
		float firstY = firstTraceHeader.y;
		boolean cdpXYFound = true;
		if(firstX == nullValue || firstY == nullValue)
		{
			float[] xy = computeAvgGatherXandY(0, 1);
			if(xy == null)
			{
				new StsMessage(frame, StsMessage.WARNING,
							   "No cdp X and Y found and Couldn't compute first trace cdp X and Y from shot/receiver records.");
				return;
			}
			cdpXYFound = false;
			firstX = xy[0];
			firstY = xy[0];
		}
		rowNumMax = lastTraceHeader.iLine;
		colNumMax = lastTraceHeader.xLine;

        float lastX, lastY;

        if(cdpXYFound)
		{
			lastX = lastTraceHeader.x;
			lastY = lastTraceHeader.y;
		}
		else
		{
			float[] xy = computeAvgGatherXandY(getNTotalTraces() - 1, -1);
			if(xy == null)
			{
				new StsMessage(frame, StsMessage.WARNING,
							   "No cdp X and Y found and Couldn't compute last cdp X and Y from shot/receiver records.");
				return;
			}
			lastX = xy[0];
			lastY = xy[0];
		}
        checkSetOriginAndAngle((double) firstX, (double) firstY, angle);
        float[] xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY((double) firstX, (double) firstY);
        xMin = xy[0];
        yMin = xy[1];
        xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY((double)lastX, (double)lastY);
        xMax = xy[0];
        yMax = xy[1];
    }

	private void isTracesPerGatherEqual(int nTracesPerGather, int nTotalTraces)
	{
		boolean isTracesPerGatherSame = true;
		// check next five gathers
		StsSEGYFormat.TraceHeader gatherLastTrace, nextGatherFirstTrace;
		if(isInline)
		{
			float currentXLine = secondTraceHeader.xLine;
			int nNextGatherFirstTrace = nTracesPerGather;
			for(int n = 0; n < 5; n++)
			{
				gatherLastTrace = getTrace(nNextGatherFirstTrace - 1);
				nextGatherFirstTrace = getTrace(nNextGatherFirstTrace);
				if(gatherLastTrace.xLine != currentXLine || nextGatherFirstTrace.xLine == currentXLine)
				{
					isTracesPerGatherSame = false;
					break;
				}
				nNextGatherFirstTrace += nTracesPerGather;
			}
		}
		else
		{
			float currentILine = secondTraceHeader.iLine;
			int nNextGatherFirstTrace = nTracesPerGather;
			for(int n = 0; n < 5; n++)
			{
				gatherLastTrace = getTrace(nNextGatherFirstTrace - 1);
				nextGatherFirstTrace = getTrace(nNextGatherFirstTrace);
				if(gatherLastTrace.iLine != currentILine || nextGatherFirstTrace.iLine == currentILine)
				{
					isTracesPerGatherSame = false;
					break;
				}
				nNextGatherFirstTrace += nTracesPerGather;
			}
		}

		if(!isTracesPerGatherSame)
		{
			new StsMessage(frame, StsMessage.WARNING,
						   "Checked first five gathers and found number of traces/gather to be variable.");
		}
		else if(nTotalTraces % nTracesPerGather != 0)
		{
			new StsMessage(frame, StsMessage.WARNING, "Number of traces per gather is variable.");
			nTracesPerGather = -1;
		}
	}

	private float[] computeAvgGatherXandY(int nStartTrace, int inc)
	{
		int nTrace = nStartTrace;
		StsSEGYFormat.TraceHeader traceHeader = getTrace(nTrace++);
		float firstIline = traceHeader.iLine;
		float firstXline = traceHeader.xLine;
		double avgX = traceHeader.x;
		double avgY = traceHeader.y;
		int nTraces = 1;
		for(int n = 0; ; n += inc)
		{
			traceHeader = getTrace(nTrace++);
			float iline = traceHeader.iLine;
			float xline = traceHeader.xLine;
			if(iline != firstIline || xline != firstXline)break;
			avgX += traceHeader.x;
			avgY += traceHeader.y;
			nTraces++;
		}
		avgX /= nTraces;
		avgY /= nTraces;
		return new float[]
			{
			(float)avgX, (float)avgY};
	}

	private int findFirstGatherTrace(byte[] tracesData, int firstTrace, int searchPointTrace)
	{
		// If we're at the first trace then no need to look further
		if (searchPointTrace == 0)
		{
			return 0;
		}
		try
		{
			StsSEGYFormat.TraceHeader firstTraceHeader, traceHeader;
			int searchPointTraceOffset = getBytesPerTrace() * (searchPointTrace - firstTrace);
			firstTraceHeader = getSegyFormat().constructTraceHeader(tracesData, searchPointTraceOffset, searchPointTrace, getIsLittleEndian());
			float firstIline = firstTraceHeader.iLine;
			float firstXline = firstTraceHeader.xLine;

			int offset = searchPointTraceOffset;
			// look backwards from the search point trace
			for (int n = searchPointTrace; n > firstTrace; n--)
			{
				offset -= getBytesPerTrace();
				traceHeader = getSegyFormat().constructTraceHeader(tracesData, offset, n - 1, getIsLittleEndian());
				float iline = traceHeader.iLine;
				float xline = traceHeader.xLine;
				if(iline != firstIline || xline != firstXline)
				{
					return n;
				}
			}
			// if we get to the start then assume that that is the first trace in the gather
			return firstTrace;
		}
		catch(Exception e)
		{
			return -1;
		}
	}
    public boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry, boolean ignoreMultiVolume, int volNum)
    {
        setNVolume(volNum, ignoreMultiVolume);
        return readWritePreStackLines(progressPanel, attributeRecords, overrideGeometry);
    }

    public boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry)
	{
		int nBlockTrace = -1;

		//this.panel = panel;
		this.attributeRecords = attributeRecords;
		nAttributes = attributeRecords.length;

        if(isCropped)
        {
            progressPanel.appendLine("PostStack3d Cropped");
        }

        if(isDepth)
		{
			byte vUnits = getSegyFormat().getZUnits();
            vertUnitsString = StsParameters.getDepthString(vUnits);
            verticalScalar = currentModel.getProject().getDepthScalar(vUnits);
		}
		else
		{
			byte vUnits = getSegyFormat().getTUnits();
            vertUnitsString = StsParameters.getTimeString(StsParameters.TIME_MSECOND);
            verticalScalar = currentModel.getProject().getTimeScalar(vUnits);
		}
		byte hUnits = getSegyFormat().getHUnits();
        this.horzUnitsString = StsParameters.getDepthString(hUnits);
        horizontalScalar = currentModel.getProject().getDepthScalar(hUnits);
        
        if(!overrideGeometry) analyzeGeometry();
		if(!initializeAttributes()) return false;
		//panel.volumeCropped.setSelected(isCropped);

		try
		{
			if(runTimer)
			{
				totalTimer = new StsTimer();
				totalTimer.start();
			}
			memoryAllocation = StsMemAllocPreStackProcess.constructor(getNTotalTraces(), getBytesPerTrace(), nSlices);
			if(memoryAllocation == null)
			{
				progressPanel.appendLine("Failed to allocate memory for SEGY volume: " + segyFilename);
				return false;
			}
			if(progressPanel.isCanceled()) return false;

			StsMessageFiles.infoMessage("Allocating " + memoryAllocation.nTotalBytes + " bytes for seismic volume processing.");

			int nTracesPerBlock = memoryAllocation.nTracesPerBlock;
			nInputBytesPerBlock = memoryAllocation.nInputBytesPerBlock;
//			nOutputSamplesPerWrite = memoryAllocation.nOutputSamplesPerWrite;
//			blockOffsetOrderBytes = memoryAllocation.blockOffsetBytes;
//			nBlocksPerWrite = memoryAllocation.nBlocksPerWrite;
			nBlockSamples = nTracesPerBlock * nSlices;
			if(debug)
			{
//				System.out.println("Input block size: " + nOutputSamplesPerWrite + ". " + nBlocksPerWrite + " input blocks per output write.");
			}
			nBlockTraces = nTracesPerBlock;
			//panel.setCurrentActionLabel("Reading SEGY volume: " + segyFilename);
			progressPanel.appendLine("Reading SEGY volume: " + segyFilename);

			if(!openFiles())return false;

			randomAccessSegyChannel = segyData.randomAccessSegyFile.getChannel();
            inputPosition = segyData.fileHeaderSize + (nVolume * segyData.nTotalTraces * segyData.bytesPerTrace);

			traceAttributes = new double[nAttributes][traceAttributesLength];

			int nBlocks = StsMath.ceiling((float)getNTotalTraces() / nTracesPerBlock);

			//panel.setCurrentActionLabel("Processing " + nBlockSamples + " samples in " + nBlocks + " blocks.");
			progressPanel.appendLine("\tProcessing volume " + (nVolume + 1) + " of " + segyData.numVolumes + " in file.");
            progressPanel.appendLine("\tProcessing " + getNTotalSamples() + " samples in " + nBlocks + " blocks.");

			Gather gather = new Gather(nSlices, getSegyFormat());
			scale = 254 / (dataMax - dataMin);
			initializeSegyIO();
			segyTraceBytes = new byte[getBytesPerTrace()];
			traceBytes = new byte[nSlices];
			int nTotalTrace = 0;
			nGatherTrace = 0;
			GatherTrace gatherTrace = null;
            checkInitializeSurveyOverride();
//			int maxNTracesTest = 30000;

//		blockLoop:
			gatherFileStartPosition = 0; // start position of buffer for next write
			gatherFileEndPosition = 0; // end position of buffer after bytes have been written

			progressPanel.initialize(nBlocks);
            progressPanel.setDescriptionAndLevel("Processing volume " + getName(), StsProgressBar.INFO);
            for(nBlock = 0; nBlock < nBlocks; nBlock++)
			{
            	Main.logUsageTimer();
				isLastBlock = (nBlock == nBlocks - 1);
				if(isLastBlock)
				{
					nBlockTraces = getNTotalTraces() - nTotalTrace;
					nInputBytesPerBlock = nBlockTraces * getBytesPerTrace();
					long remaining = getSegyFileSize() - inputPosition;
					System.out.println("Input file remaining bytes " + remaining + " need bytes " + nInputBytesPerBlock);

				}
				// We may have traces from previous gather which haven't been written out so
				// buffer size needs to be size of this block plus what hasn't been written yet.
				// So just bump the end position by the block size and compute the total size.
				gatherFileStartPosition = nGatherTracesWritten * nSlices;
				gatherFileEndPosition = (nLineTraces + nBlockTraces) * nSlices;
				long gatherBufferSize = gatherFileEndPosition - gatherFileStartPosition;
				gatherFloatBuffer.clean();
				if(!gatherFloatBuffer.map(gatherFileStartPosition, gatherBufferSize))
				{
					closeGatherFile();
					deletePartialFiles();
					return false;
				}
				if(debugLine)System.out.println("Mapping gather output floatBuffer start position: " + gatherFileStartPosition + " size: " + gatherBufferSize);
//                runTimer = true;
				if(runTimer)timer.start();
				mapInputBuffer();
//				mapOutputBuffers();
				for(nBlockTrace = 0; nBlockTrace < nBlockTraces; nBlockTrace++, nTotalTrace++, nLineTraces++, nGatherTrace++)
				{
//					if(nTotalTrace >= maxNTracesTest) break blockLoop;
//					System.out.println("Block " + nBlock + " trace " + nTrace);
					inputBuffer.get(segyTraceBytes);
					GatherTrace prevGatherTrace = gatherTrace;
					gatherTrace = gather.addGatherTrace(segyTraceBytes);
					gather.checkTrace(gatherTrace, prevGatherTrace);

					if(nLineTraces == 0) // only true the first time thru; when a line is finished below, the next gather will have been found and nLineTraces will be 1
					{
						initializeLine();
						gather.startLine(overrideGeometry);
					}
					if(gather.gatherChanged) // end this gather
					{
						gather.outputGather(nGatherTrace, true);
						setLastTraceInGather(nLineTraces - 1);

                        // Missing Gather Processing
                        if(!gather.lineChanged && (gatherTrace.xline-prevGatherTrace.xline) > colNumInc)
                        {
                            float tempXline = prevGatherTrace.xline + colNumInc;
                            while(tempXline < gatherTrace.xline)
                            {
                                setLastTraceInGather(nLineTraces - 1);
                                tempXline = tempXline + colNumInc;
                            }
                            gather.gatherChanged = true;
                        }

                    }
					if(gather.lineChanged) // end this line
					{
						gather.endLine(prevGatherTrace, overrideGeometry);
                        writeLineFiles();
                        progressPanel.appendLine("\tSequential line " + nLine + " output " + nLineTraces + " traces in " + nCols + " gathers.");
					}
					if(gather.gatherChanged) // start next gather
					{
						gather.startNextGather(nGatherTrace);
					}
					if(gather.lineChanged) // start next line
					{
						nLine++;
						openEmptyGatherFile();

//						mapOutputBuffers();
						gather.lineChanged = false;
						initializeLine();

						gatherFileStartPosition = 0;
						long gatherFileSize = nBlockTraces * nSlices;
						nGatherTracesWritten = 0;
						gatherFileEndPosition = gatherFileStartPosition + gatherFileSize;
						if(!gatherFloatBuffer.map(0, gatherFileSize)) return false;
						if(debugLine)System.out.println("Mapping gather output floatBuffer. position: 0" + " size: " + nBlockTraces * nSlices);

						gather.startLine(overrideGeometry);
					}
				}
				clearInputBuffer();
//				clearOutputBuffers();
				if(runTimer && (nBlock+1)%10 == 0) totalTimer.stopPrint("Time to output 10 blocks thru block " + nBlock);
                if(progressPanel.isCanceled())
				{
					progressPanel.appendLine("\tProcess interrupted. Deleting incomplete lines.");

					// delete current files to prvent them being left in an incoherent state
					closeGatherFile();
					deletePartialFiles();
					return false;
				}
                progressPanel.setValue(nBlock);
            }

			gather.outputLastGather();
			gather.endLine(gatherTrace, overrideGeometry);
            setLastTraceInGather(nLineTraces-1);
            writeLineFiles();
            if(memoryAllocation != null)
            {
                System.out.println("Freeing memory.....");
//                memoryAllocation.freeBuffers();
                memoryAllocation = null;
            }
			progressPanel.appendLine("\tSEGY PostStack3d Read");
			progressPanel.appendLine("\tGather PostStack3d Written");
			progressPanel.appendLine("\tVolume Header File Written");
//            progressPanel.appendLine("\tFinished in: " + totalTimer.getElapsedTime()/1000 + " seconds\n");

			//panel.cubeRead.setSelected(true);
			//panel.gatherCubeWritten.setSelected(true);
//rjcxx			calculateHistogram();
			//panel.headerWritten.setSelected(true);
			//panel.setCurrentActionLabel("");
			if(runTimer) totalTimer.stopPrint("Time to output all files:");
//            runTimer = false;
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackSegyLine3d.readWritePreStackLines() failed.", e, StsException.WARNING);
			return false;
		}
	}

    protected boolean writeLineFiles()
    {
        closeGatherFile();
        renameGathersFile();
        trimGathersArray();
        writeAttributesFile();
        writeHeaderFile();
        return true;
    }

    public void renameGathersFile()
    {
        File file = new File(stsDirectory, gatherFilename);
        String newName = createFilename(gatherFormat);
        if(gatherFilename.equals(newName)) return;
        gatherFilename = newName;
        File newFile = new File(stsDirectory, newName);
        file.renameTo(newFile);
    }

    public boolean initializeAttributes()
	{
		int nAttributes = attributeRecords.length;
		attributeNames = new String[nAttributes];
		attributeMinValues = new double[nAttributes];
		attributeMaxValues = new double[nAttributes];
		for(int n = 0; n < nAttributes; n++)
		{
			StsSEGYFormatRec record = attributeRecords[n];
			String userName = record.getUserName();
			attributeNames[n] = userName;
			attributeMinValues[n] = StsParameters.largeDouble;
			attributeMaxValues[n] = -StsParameters.largeDouble;
//			if (userName.equals("CDP"))
//				cdpAttributeIndex = n;
			if(userName.equals("ILINE_NO"))
				ilineAttributeIndex = n;
			else if(userName.equals("XLINE_NO"))
				xlineAttributeIndex = n;
			else if(userName.equals(StsSEGYFormat.CDP_X))
				cdpxAttributeIndex = n;
			else if(userName.equals(StsSEGYFormat.CDP_Y))
				cdpyAttributeIndex = n;
			else if(userName.equals("OFFSET"))
				offsetAttributeIndex = n;
		}
		boolean attributesOk = true;
		/*
		 if(cdpAttributeIndex == -1)
		 {
		  attributesOk = false;
		  attributeMissingMessage("trace CDP number");
		 }
		 */
		if(ilineAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace inline number");
		}
		if(xlineAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace crossline number");
		}
		if(cdpxAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace CDP X value");
		}
		if(cdpyAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace CDP Y value");
		}
		if(offsetAttributeIndex == -1)
		{
			attributesOk = false;
			attributeMissingMessage("trace offset value");
		}
		return attributesOk;
	}

	private void attributeMissingMessage(String name)
	{
		new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + name + " attribute");
	}

	public void initializeLine()
	{
		nLastTraceInGathers = new int[10000];
		nCols = 0;
		nOffsetsMax = 0;
		nLineTraces = 0;
		traceOffsetMin = StsParameters.largeFloat;
		traceOffsetMax = -StsParameters.largeFloat;
	}

    private void checkInitializeSurveyOverride()
    {
		if(this.overrideGeometryBox == null) return;
        xOrigin = overrideGeometryBox.xOrigin;
        yOrigin = overrideGeometryBox.yOrigin;
        rowNumInc = overrideGeometryBox.rowNumInc;
        colNumInc = overrideGeometryBox.colNumInc;
        originRowNum = overrideGeometryBox.rowNumMin;
        originColNum = overrideGeometryBox.colNumMin;
        angle = overrideGeometryBox.angle;
        xInc = overrideGeometryBox.xInc;
        yInc = overrideGeometryBox.yInc;
        xMin = overrideGeometryBox.xMin;
        yMin = overrideGeometryBox.yMin;
        zMin = overrideGeometryBox.zMin;
        zInc = overrideGeometryBox.zInc;
        nSlices = overrideGeometryBox.nSlices;
        zMax = zMin + zInc*(nSlices-1);
        isXLineCCW = overrideGeometryBox.isXLineCCW;
    }

	protected boolean areTracesInSameGather(TraceHeader firstTraceHeader, TraceHeader secondTraceHeader)
	{
		return (firstTraceHeader.iLine == secondTraceHeader.iLine && firstTraceHeader.xLine == secondTraceHeader.xLine);
	}

	public int compareTo(Object o)
	{
		return 0;
	}

	public class Gather
	{
		int nSamples;
		StsSEGYFormat segyFormat = null;
		GatherTrace[] gatherTraces = new GatherTrace[100];
		public boolean lineChanged = false;
		public boolean gatherChanged = false;
		public float iline;
		public float xline;

		public Gather(int nSamples, StsSEGYFormat segyFormat)
		{
			this.segyFormat = segyFormat;
			this.nSamples = nSamples;
			//           initializeAttributes();
		}

		public GatherTrace addGatherTrace(byte[] segyTraceBytes)
		{
			int length = gatherTraces.length;
			if(nGatherTrace > length - 1)
			{
				GatherTrace[] newGatherTraces = new GatherTrace[length + 100];
				System.arraycopy(gatherTraces, 0, newGatherTraces, 0, length);
				gatherTraces = newGatherTraces;
			}

			GatherTrace gatherTrace = gatherTraces[nGatherTrace];
			if(gatherTrace == null)
			{
				gatherTrace = new GatherTrace(nAttributes, nSamples);
				gatherTraces[nGatherTrace] = gatherTrace;
			}

			gatherTrace.processTraceAttributes(segyTraceBytes, attributeRecords, this);
			gatherTrace.processTrace(segyTraceBytes, nSamples);
			if(nGatherTrace == 0)
			{
				iline = gatherTrace.iline;
				xline = gatherTrace.xline;
//				cdp = gatherTrace.cdp;
			}
			return gatherTrace;
		}

		void startNextGather(int nNextTrace)
		{
			// move last trace which is first trace of next line to position 0; swap traces so there are no copies in array
			GatherTrace gatherTrace0 = gatherTraces[0];
			gatherTraces[0] = gatherTraces[nNextTrace]; // last trace is first of next line, so move there
			gatherTraces[nNextTrace] = gatherTrace0;
			nGatherTrace = 0;
			iline = gatherTraces[0].iline;
			xline = gatherTraces[0].xline;
//			cdp = gatherTraces[0].cdp;
			gatherChanged = false;
		}

		public void startLine(boolean override)
		{
			double[] attributes = gatherTraces[0].attributes;
			iline = gatherTraces[0].iline;
			xline = gatherTraces[0].xline;
            rowNumMin = iline;
            colNumMin = xline;
            rowNumMax = iline;
            colNumMax = xline;
            if(override)
            {
                if(!overrideGeometryBox.isXLineCCW)
                    yMin = overrideGeometryBox.yMin - yInc*(iline - originRowNum)/rowNumInc;
                else
                    yMin = overrideGeometryBox.yMin + yInc*(iline - originRowNum)/rowNumInc;
                xMin = overrideGeometryBox.xMin + xInc*(xline - originColNum)/colNumInc;
                double[] absXY = getAbsoluteXY(xMin, yMin);
                firstX = (float)absXY[0];
                firstY = (float)absXY[1];
			}
			else
			{
				firstX = (float)attributes[cdpxAttributeIndex];
				firstY = (float)attributes[cdpyAttributeIndex];
				checkSetOriginAndAngle((double) firstX, (double) firstY, angle);
				float[] xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY((double) firstX, (double) firstY);
				xMin = xy[0];
				yMin = xy[1];
                System.out.println("Line Mins: Y=" + yMin + " X=" + xMin + " Iline=" + iline + " xLine=" + xline);
			}
		}

		public void checkTrace(GatherTrace gatherTrace, GatherTrace prevGatherTrace)
		{
			if(prevGatherTrace == null)return;
			if(isInline)
			{
				if(gatherTrace.iline != rowNumMin)
					lineChanged = true;
				if(gatherTrace.xline != xline)
					gatherChanged = true;
			}
			else
			{
				if(gatherTrace.xline != colNumMin)
					lineChanged = true;
				if(gatherTrace.iline != iline)
					gatherChanged = true;
			}
			if(!gatherChanged && lineChanged)
				StsException.systemError("Line changed but gather not changed.");
		}

		public void endLine(GatherTrace gatherTrace, boolean override)
		{
			double[] attributes = gatherTrace.attributes;
			iline = gatherTrace.iline;
			xline = gatherTrace.xline;
            rowNumMax = iline;
            colNumMax = xline;
            if(override)
			{
                if(!overrideGeometryBox.isXLineCCW)
                    yMax = overrideGeometryBox.yMin - yInc*(iline - originRowNum)/rowNumInc;
                else
                    yMax= overrideGeometryBox.yMin + yInc*(iline - originRowNum)/rowNumInc;
                xMax = overrideGeometryBox.xMin + xInc*(xline - originColNum)/colNumInc;
				double[] absXY = getAbsoluteXY(xMin, yMin);
                lastX = (float)absXY[0];
                firstY = (float)absXY[1];
			}
			else
			{
				lastX = (float)attributes[cdpxAttributeIndex];
				lastY = (float)attributes[cdpyAttributeIndex];
				float[] xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY((double)lastX, (double)lastY);
				xMax = xy[0];
				yMax = xy[1];
//                System.out.println("Line Maxs: Y=" + yMax + " X=" + xMax);
			}

		}

		public void outputLastGather()
		{
			outputGather(nGatherTrace, false);
//			outputGather(nGatherTrace + 1, false);
		}

		// called if we've changed gathers and we only put out the previous gather (not including this new trace)
		public void outputGather(int nGatherTraces, boolean checkGather)
		{
            outputGather(nGatherTraces, checkGather, false);
		}

        public void outputGather(int nGatherTraces, boolean checkGather, boolean deadTrace)
        {
			if(nGatherTraces == 0) return;
			nOffsetsMax = Math.max(nOffsetsMax, nGatherTraces);
            if(debugGather)
			{
				System.out.println("Line " + nLine + " gather " + nCols + " sort and output " + nGatherTraces + " traces.");
				System.out.println("         first " + gatherTraces[0].iline + " " + gatherTraces[0].xline + " last " +
								   gatherTraces[nGatherTraces - 1].iline + " " + gatherTraces[nGatherTraces - 1].xline);
                System.out.println("nGathers=" + nCols + " nGatherTraces=" + nGatherTraces + " nGatherTracesWritten=" + nGatherTracesWritten);

            }
			Arrays.sort(gatherTraces, 0, nGatherTraces);

			if(nLineTraces >= traceAttributesLength - 1) increaseTraceAttributesArray();
			for(int n = 0; n < nGatherTraces; n++)
			{
                gatherTraces[n].output(deadTrace);
				gatherTraces[n].saveAttributes();
			}
			nGatherTracesWritten += nGatherTraces;
			if(debugGather)
			{
//				System.out.println("Current floatBuffer position: " + nGatherTracesWritten*nCroppedSlices + " Should be " + nLineTraces * nCroppedSlices);
//				System.out.println("Free memory is "+Runtime.getRuntime().freeMemory()+" total "+Runtime.getRuntime().totalMemory());
			}
		}
	}

	public class GatherTrace implements Comparable
	{
		double[] attributes;
		float[] values, deadTrace;
		float traceOffset;
		float iline, xline;
//		int cdp;
//		double[] values;
//		float offset;

		public GatherTrace(int nAttributes, int nSamples)
		{
			attributes = new double[nAttributes];
			values = new float[nSamples];
            deadTrace = new float[nSamples];
            for(int i=0; i<nSamples; i++)
                deadTrace[i] = 0.0f;
//			values = new double[nSamples];
		}

		private boolean processTraceAttributes(byte[] traceBytes, StsSEGYFormatRec[] attributeRecords, Gather gather)
		{
			int nAttributes = attributeRecords.length;
            float edScale = getSegyFormat().getFloatElevationScale(traceBytes, getIsLittleEndian());
            float xyScale = getSegyFormat().getFloatCoordinateScale(traceBytes, getIsLittleEndian());

			for(int n = 0; n < nAttributes; n++)
			{
				StsSEGYFormatRec record = attributeRecords[n];
				int loc = record.getLoc();
				int format = record.getFormat();

				float appliedScale = 1.0f;
				if(record.getApplyScalar().equals("CO-SCAL"))
					appliedScale = xyScale;
                else if(record.getApplyScalar().equals("ED-SCAL"))
                    appliedScale = edScale;

				switch(format)
				{
					case StsSEGYFormat.IBMFLT: // IBM Float
						attributes[n] = (float)StsMath.convertIBMFloatBytes(traceBytes, loc, getIsLittleEndian()) * appliedScale;
						break;
					case StsSEGYFormat.IEEEFLT: // IEEE Float
						attributes[n] = (float)Float.intBitsToFloat(StsMath.convertIntBytes(traceBytes, loc, getIsLittleEndian())) * appliedScale;
						break;
					case StsSEGYFormat.INT4: // Integer 4
						attributes[n] = (float)StsMath.convertIntBytes(traceBytes, loc, getIsLittleEndian()) * appliedScale;
						break;
					case StsSEGYFormat.INT2: // Integer 2
						attributes[n] = (float)StsMath.convertBytesToShort(traceBytes, loc, getIsLittleEndian()) * appliedScale;
						break;
					default:
						;
				}
			}
			traceOffset = (float)attributes[offsetAttributeIndex];
			iline = (float)attributes[ilineAttributeIndex];
			xline = (float)attributes[xlineAttributeIndex];
//			cdp = (int)attributes[cdpAttributeIndex];

			traceOffsetMin = Math.min(traceOffsetMin, traceOffset);
			traceOffsetMax = Math.max(traceOffsetMax, traceOffset);
			/*
					 for(int n = 0; n < nAttributes; n++)
			   {
				traceAttributes[n][nLineTrace] = attributes[n];
			   }
			 */
			return true;
		}

		void saveAttributes()
		{
			try
			{
				for(int n = 0; n < nAttributes; n++)
				{
					double attribute = attributes[n];
					traceAttributes[n][nAttributeTrace] = attribute;
					attributeMinValues[n] = Math.min(attributeMinValues[n], attribute);
					attributeMaxValues[n] = Math.max(attributeMaxValues[n], attribute);
				}
				nAttributeTrace++;
			}
			catch(Exception e)
			{
				StsException.systemError("StsPreStackSegyLine3d.saveAttributes() failed.");
			}
		}

		void processTrace(byte[] segyTraceBytes, int nSamples)
		{
			segyIO.processTrace(segyTraceBytes, nSamples, nLineTraces, values);
//			segyIO.processTraceDouble(segyTraceBytes, nSamples, nTrace, values);

		}

		public int compareTo(Object other)
		{
			float otherOffset = ((GatherTrace)other).traceOffset;
			if(traceOffset > otherOffset)return 1;
			if(traceOffset < otherOffset)return -1;
			return 0;
		}

        void output()
		{
            output(false);
        }
		void output(boolean dead)
		{
			try
			{
				if(!dead)
                   gatherFloatBuffer.put(values);
                else
                   gatherFloatBuffer.put(deadTrace);
//			    for(int n = 0; n < nCroppedSlices; n++)
//					gatherDOS.writeFloat(values[n]);
			}
			catch(Exception e)
			{
				StsException.outputException("StsPreStackSegyLine3d.GatherTrace.output() failed.", e, StsException.WARNING);
			}
		}
	}

}