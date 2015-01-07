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
public  class StsPreStackSegyLine2d extends StsPreStackSegyLine
{
	/** indexes for key attribute records */
	transient int cdpAttributeIndex = -1;

    public String getGroupname()
    {
        return group2dPrestack;
    }

    public StsPreStackSegyLine2d()
	{
		super();
	}

	public StsPreStackSegyLine2d(boolean persistent)
	{
		super(persistent);
	}

	protected StsPreStackSegyLine2d(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws FileNotFoundException, StsException
	{
		super(seismicWizard, file, stsDirectory, frame, segyFormat);
	}

	static public StsPreStackSegyLine2d constructor(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat)
	{
		try
		{
			return new StsPreStackSegyLine2d(seismicWizard, file, stsDirectory, frame, segyFormat);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSegyline2d.constructTraceAnalyzer(file, stsDirectory) failed.", e,
										 StsException.WARNING);
			return null;
		}
	}

    /** computes line geometry from trace CDP i X, and Y values */
	public boolean analyzeGeometry()
	{
		try
		{
            originSet = false;
			if(segyData.nTotalTraces == 0) return false;

			firstTraceHeader = getTrace(0);
			if(firstTraceHeader == null)return false;


			lastTraceHeader = getTrace(getNTotalTraces() - 1);
			if(lastTraceHeader == null) return false;

			xOrigin = firstTraceHeader.x;
			yOrigin = firstTraceHeader.y;

            rowNumMin = firstTraceHeader.iLine;
			colNumMin = firstTraceHeader.xLine;
			rowNumMax = lastTraceHeader.iLine;
			colNumMax = lastTraceHeader.xLine;

//			rowNumMin = nullValue;
//			colNumMin = nullValue;
//			rowNumMax = nullValue;
//			colNumMax = nullValue;
			rowNumInc = 0;
			colNumInc = 0;

            secondTraceHeader = findFirstTraceSecondGather();
			if(secondTraceHeader == null)return false;
			return true;
		}
		catch(Exception e)
		{
			StsMessageFiles.errorMessage("Couldn't analyze geometry: select a different number of bytes per sample.");
			return false;
		}
	}

	public double getEndX()
	{
		if (lastTraceHeader == null)
			return 0d;
		return lastTraceHeader.x;
	}

	public double getEndY()
	{
		if (lastTraceHeader == null)
			return 0d;
		return lastTraceHeader.y;
	}

	private StsSEGYFormat.TraceHeader findFirstTraceSecondGather()
	{
		secondTraceHeader = null;

		int cdpMin = firstTraceHeader.cdp;
		nTracesPerGather = 1;
		while(nTracesPerGather < 1000)
		{
			secondTraceHeader = getTrace(nTracesPerGather++);
			if(secondTraceHeader.cdp != cdpMin)
			{
				break;
			}
			secondTraceHeader = null;
		}
		if(secondTraceHeader == null)
		{
			StsMessageFiles.errorMessage("Unable to find next gather.");
			return null;
		}
		return secondTraceHeader;
	}

	protected boolean areTracesInSameGather(TraceHeader firstTraceHeader, TraceHeader secondTraceHeader)
	{
		return (firstTraceHeader.cdp == secondTraceHeader.cdp);
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

    public boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry, boolean ignoreMultiVolume, int volumeNum)
    {
        setNVolume(volumeNum, ignoreMultiVolume);
        return readWritePreStackLines(progressPanel, attributeRecords, overrideGeometry);
    }

    public boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry)
    {
		int nBlockTrace = -1;

		//this.panel = panel;
		this.attributeRecords = attributeRecords;
		nAttributes = attributeRecords.length;

		if(isCropped) progressPanel.appendLine("PostStack3d Cropped");

		if(isDepth)
		{
			byte vUnits = getSegyFormat().getZUnits();
			verticalScalar = currentModel.getProject().getDepthScalar(vUnits);
		}
		else
		{
			byte vUnits = getSegyFormat().getTUnits();
			verticalScalar = currentModel.getProject().getTimeScalar(vUnits);
		}
		byte hUnits = getSegyFormat().getHUnits();
		horizontalScalar = currentModel.getProject().getDepthScalar(hUnits);
		if(!overrideGeometry)
           analyzeGeometry();
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
				progressPanel.appendLine("\tFailed to allocate memory for SEGY line: " + segyFilename);
				return false;
			}
			if(progressPanel.isCanceled()) return false;

			StsMessageFiles.infoMessage("Allocating " + memoryAllocation.nTotalBytes +
										" bytes for seismic line processing.");

			int nTracesPerBlock = memoryAllocation.nTracesPerBlock;
			nInputBytesPerBlock = memoryAllocation.nInputBytesPerBlock;
			nBlockSamples = nTracesPerBlock * nSlices;
			if(debug)
			{
//				System.out.println("Input block size: " + nOutputSamplesPerWrite + ". " + nBlocksPerWrite +
//								   " input blocks per output write.");
			}
			nBlockTraces = nTracesPerBlock;

			progressPanel.appendLine("\tReading SEGY line: " + segyFilename);

			if(!openEmptyGatherFile())return false;

			randomAccessSegyChannel = segyData.randomAccessSegyFile.getChannel();
			inputPosition = segyData.fileHeaderSize;

			traceAttributes = new double[nAttributes][traceAttributesLength];

			int nBlocks = StsMath.ceiling((float)getNTotalTraces() / nTracesPerBlock);

			//panel.setCurrentActionLabel("Processing " + nBlockSamples + " samples in " + nBlocks + " blocks.");
			progressPanel.appendLine("\tProcessing " + getNTotalSamples() + " samples in " + nBlocks + " blocks.");

			Gather gather = new Gather(nSlices, getSegyFormat());
			scale = 254 / (dataMax - dataMin);
			initializeSegyIO();
			segyTraceBytes = new byte[getBytesPerTrace()];
			traceBytes = new byte[nSlices];
			int nTotalTrace = 0;
			nGatherTrace = 0;
			GatherTrace gatherTrace = null;
//			int maxNTracesTest = 30000;

//		blockLoop:
			gatherFileStartPosition = 0; // start position of buffer for next write
			gatherFileEndPosition = 0; // end position of buffer after bytes have been written

			progressPanel.initialize(nBlocks);
            progressPanel.setDescriptionAndLevel("Processing 2d line " + getName(), StsProgressBar.INFO);
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
				if(!gatherFloatBuffer.map(gatherFileStartPosition, gatherBufferSize, progressPanel)) return false;
				if(debugLine)System.out.println("Mapping gather output floatBuffer start position: " + gatherFileStartPosition + " size: " + gatherBufferSize);
				if(runTimer)timer.start();
				mapInputBuffer();
				for(nBlockTrace = 0; nBlockTrace < nBlockTraces; nBlockTrace++, nTotalTrace++, nLineTraces++, nGatherTrace++)
				{
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
//						setGatherXY(prevGatherTrace);
						setLastTraceInGather(nLineTraces - 1);

                        // Missing Gather Processing ---- NOT TEST - DID NOT HAVE DATASET
                        if((gatherTrace.cdp-prevGatherTrace.cdp) > 1)
                        {
                            int tempCDP = prevGatherTrace.cdp + 1;
                            while(tempCDP < gatherTrace.cdp)
                            {
                                setLastTraceInGather(nLineTraces - 1);
                                tempCDP = tempCDP + 1;
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
						openFiles();
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
				if(runTimer)timer.stopPrint("Time to output block " + nBlock);
				if (progressPanel.isCanceled())
				{
					progressPanel.appendLine("\tProcess interrupted. Deleting incomplete lines.");
					closeGatherFile();
					deletePartialFiles();
					return false;
				}
                progressPanel.setValue(nBlock+1);
            }

			gather.outputLastGather();
			gather.endLine(gatherTrace, overrideGeometry);
            setLastTraceInGather(nLineTraces-1);
            writeLineFiles();

			progressPanel.appendLine("\tSEGY 2d Line Read");
			progressPanel.appendLine("\t2d Line Gathers File Written");
			progressPanel.appendLine("\t2d Line Header File Written");
            progressPanel.appendLine("\t2d Line Attributes File Written");
			if(runTimer) totalTimer.stopPrint("Time to output all files:");
//            runTimer = false;
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackSegyLine3d.readWritePreStackLines() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("Exception thrown: " + e.getMessage(), StsProgressBar.ERROR);
            progressPanel.finished();
            return false;
		}
	}

    public void renameGathersFile()
    {
    }
    private boolean initializeAttributes()
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
			if (userName.equals(StsSEGYFormat.CDP))
				cdpAttributeIndex = n;
			if(userName.equals("ILINE_NO"))
				ilineAttributeIndex = n;
			else if(userName.equals("XLINE_NO"))
				xlineAttributeIndex = n;
			else if(userName.equals(StsSEGYFormat.CDP_X))
				cdpxAttributeIndex = n;
			else if(userName.equals(StsSEGYFormat.CDP_Y))
				cdpyAttributeIndex = n;
            else if(userName.equals(StsSEGYFormat.SHT_X))
                shtxAttributeIndex = n;
            else if(userName.equals(StsSEGYFormat.SHT_Y))
                shtyAttributeIndex = n;
            else if(userName.equals(StsSEGYFormat.REC_X))
                recxAttributeIndex = n;
            else if(userName.equals(StsSEGYFormat.REC_Y))
                recyAttributeIndex = n;
			else if(userName.equals("OFFSET"))
				offsetAttributeIndex = n;
		}
		boolean attributesOk = true;

		 if(cdpAttributeIndex == -1)
		 {
		  attributesOk = false;
		  attributeMissingMessage("trace CDP number");
		 }
		/*
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
	    */
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
/*
	private void setGatherXY(GatherTrace gatherTrace)
	{
		cdpX[nGathers] = gatherTrace.cdpX;
		// jbw test put a kink in it
		//if (nGathers > 500)
		//cdpY[nGathers] = gatherTrace.cdpY + (33. * (nGathers - 500));
		//else
		cdpY[nGathers] = gatherTrace.cdpY ;
		cdp[nGathers]= gatherTrace.cdp;
	    //cdpY[nGathers] -=500.;
    }
 */
	class Gather
	{
		int nSamples;
		StsSEGYFormat segyFormat = null;
		GatherTrace[] gatherTraces = new GatherTrace[100];
		boolean lineChanged = false;
		boolean gatherChanged = false;
		//float iline;
		//float xline;
		int cdp;

		Gather(int nSamples, StsSEGYFormat segyFormat)
		{
			this.segyFormat = segyFormat;
			this.nSamples = nSamples;
			//           initializeAttributes();
		}

		GatherTrace addGatherTrace(byte[] segyTraceBytes)
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
				//iline = gatherTrace.iline;
				//xline = gatherTrace.xline;
				cdp = gatherTrace.cdp;
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
			//iline = gatherTraces[0].iline;
			//xline = gatherTraces[0].xline;
			cdp = gatherTraces[0].cdp;
			gatherChanged = false;
		}

		void startLine(boolean override)
		{
            if(overrideGeometryBox == null)
            {
                double[] attributes = gatherTraces[0].attributes;
                xOrigin = attributes[cdpxAttributeIndex];
                yOrigin = attributes[cdpyAttributeIndex];
            }
        }

		void checkTrace(GatherTrace gatherTrace, GatherTrace prevGatherTrace)
		{
			if(prevGatherTrace == null)return;
			if (gatherTrace.cdp != cdp)
				gatherChanged = true;
			return;
		}

		void endLine(GatherTrace gatherTrace, boolean override)
		{
            if(overrideGeometryBox != null) return;

            double[] attributes = gatherTrace.attributes;
			double x = attributes[cdpxAttributeIndex];
		    double y = attributes[cdpyAttributeIndex];

            // Construct our local coordinate boundingBox around the data using the first and last points.
            // If the line bends, we possibly don't have a box around all the data, which only means the 3d
            // cursor won't cover the full extent.  Fix if necessary by getting the boundingBox for all the data.
            // Given the first and last XYs, the absolute XY origin remains at the first XY.
            // If the last X is less than the first X, then the local X will run from -dX to 0.  Same logic for Y.

            double dx = x - xOrigin;
            double dy = y - yOrigin;

            if(dx < 0.0)
            {
                xMin = (float)dx;
                xMax = 0.0f;
            }
            else
            {
                xMin = 0.0f;
                xMax = (float)dx;
            }

            if(dy < 0.0)
            {
                yMin = (float)dy;
                yMax = 0.0f;
            }
            else
            {
                yMin = 0.0f;
                yMax = (float)dy;
            }
		}

		void outputLastGather()
		{
			outputGather(nGatherTrace, false);
		}

		// called if we've changed gathers and we only put out the previous gather (not including this new trace)
        void outputGather(int nGatherTraces, boolean checkGather)
        {
            outputGather(nGatherTraces, checkGather, false);
        }

		void outputGather(int nGatherTraces, boolean checkGather, boolean deadTrace)
		{
			if(nGatherTraces == 0)return;
			nOffsetsMax = Math.max(nOffsetsMax, nGatherTraces);
			if(debugGather)
			{
				System.out.println("Line " + nLine + " gather " + nCols + " sort and output " + nGatherTraces + " traces.");
				System.out.println("         first " + gatherTraces[0].cdp + " " + gatherTraces[0].cdp + " last " +
								   gatherTraces[nGatherTraces - 1].cdp + " " + gatherTraces[nGatherTraces - 1].cdp);
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

	class GatherTrace implements Comparable
	{
		double[] attributes;
		float[] values, deadTrace;
		float traceOffset;
		//float iline, xline;
		int cdp;
		double cdpX;
		double cdpY;
//		double[] values;
//		float offset;

		GatherTrace(int nAttributes, int nSamples)
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
				if(record.getApplyScalar().equals("CO-SCAL"))  // jbw changed from !
					appliedScale = xyScale;
				else if(record.getApplyScalar().equals("ED-SCAL")) // jbw
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
			//iline = (float)attributes[ilineAttributeIndex];
			//xline = (float)attributes[xlineAttributeIndex];
			cdp = (int)attributes[cdpAttributeIndex];
			cdpX = attributes[cdpxAttributeIndex];
			cdpY = attributes[cdpyAttributeIndex];
            if((cdpX == 0.0) && (cdpY == 0.0))
            {
                cdpX = calculateCdpX();
                cdpY = calculateCdpY();
            }
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

        private double calculateCdpX()
        {
            if((shtxAttributeIndex == -1) || (recxAttributeIndex == -1))
                return 0.0f;
            else
            {
                attributes[cdpxAttributeIndex] = attributes[shtxAttributeIndex] + attributes[recxAttributeIndex];
                return attributes[cdpxAttributeIndex];
            }
        }
        private double calculateCdpY()
        {
            if((shtyAttributeIndex == -1) || (recyAttributeIndex == -1))
                return 0.0f;
            else
            {
                attributes[cdpyAttributeIndex] = attributes[shtyAttributeIndex] + attributes[recyAttributeIndex];
                return attributes[cdpyAttributeIndex];
            }
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
