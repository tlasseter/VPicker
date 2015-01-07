package com.Sts.Types;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;

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
public  class StsPreStackMicroseismic extends StsPreStackSegyLine3d
{
    public StsPreStackMicroseismic()
	{
		super();
	}

	public StsPreStackMicroseismic(boolean persistent)
	{
		super(persistent);
	}

	protected StsPreStackMicroseismic(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat) throws FileNotFoundException, StsException
	{
		super(seismicWizard, file, stsDirectory, frame, segyFormat);
	}

	static public StsPreStackMicroseismic constructor(StsSeismicWizard seismicWizard, StsFile file, String stsDirectory, Frame frame, StsSEGYFormat segyFormat)
	{
		try
		{
			return new StsPreStackMicroseismic(seismicWizard, file, stsDirectory, frame, segyFormat);
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackMicroseismic.constructor(seismicWizard, file, stsDirectory, frame, segyFormat) failed.", e,
										 StsException.WARNING);
			return null;
		}
	}
    /** computes volume geometry from trace CDP ilines, crossline, X, and Y values */
	public boolean analyzeGeometry()
	{
		try
		{
            originSet = false;

			// Get first trace and first trace in next gather to determine if this is line or xline
			firstTraceHeader = getTrace(0);
			if(firstTraceHeader == null)return false;

			lastTraceHeader = getTrace(getNTotalTraces() - 1);
			if(lastTraceHeader == null)return false;

			xOrigin = firstTraceHeader.x;
			yOrigin = firstTraceHeader.y;

			rowNumMin = 1;
			colNumMin = firstTraceHeader.xLine;
			rowNumMax = 1;
			colNumMax = lastTraceHeader.xLine;

			findLastTraceFirstLine();
			findFirstTraceLastLine();

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
	
    // TODO this method is very similar to method in StsPreStackSegyVsp, so with a little
    // TODO effort, it could be moved up to StsPreStackSegyLine and used by both with some common working methods in each subClass
    public boolean readWritePreStackLines(StsProgressPanel progressPanel, StsSEGYFormatRec[] attributeRecords, boolean overrideGeometry)
	{
		int nBlockTrace = -1;

		//this.panel = panel;
		this.attributeRecords = attributeRecords;
		nAttributes = attributeRecords.length;

		byte vUnits = getSegyFormat().getTUnits();
		verticalScalar = currentModel.getProject().getTimeScalar(vUnits);
			
		byte hUnits = getSegyFormat().getHUnits();
		horizontalScalar = currentModel.getProject().getDepthScalar(hUnits);
		if(!overrideGeometry) analyzeGeometry();
		if(!initializeAttributes()) return false;
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
			nBlockSamples = nTracesPerBlock * nSlices;
			nBlockTraces = nTracesPerBlock;
			progressPanel.appendLine("Reading SEGY gather: " + segyFilename);

			if(!openFiles())return false;

			randomAccessSegyChannel = segyData.randomAccessSegyFile.getChannel();
            inputPosition = segyData.fileHeaderSize + (nVolume * segyData.nTotalTraces * segyData.bytesPerTrace);

			traceAttributes = new double[nAttributes][traceAttributesLength];

			int nBlocks = StsMath.ceiling((float)getNTotalTraces() / nTracesPerBlock);

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

			gatherFileStartPosition = 0; // start position of buffer for next write
			gatherFileEndPosition = 0; // end position of buffer after bytes have been written

			progressPanel.initialize(nBlocks);
            progressPanel.setDescriptionAndLevel("Processing gather " + getName(), StsProgressBar.INFO);
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
					if(gather.lineChanged) // end this line
					{
						gather.outputGather(nGatherTrace, true);
						setLastTraceInGather(nLineTraces - 1);						
						gather.endLine(prevGatherTrace, overrideGeometry);
                        writeLineFiles();
                        progressPanel.appendLine("\tSequential line " + nLine + " output " + nLineTraces + " traces in " + nCols + " gathers.");
					}
				}
				clearInputBuffer();
				if(runTimer && (nBlock+1)%10 == 0) totalTimer.stopPrint("Time to output 10 blocks thru block " + nBlock);
                if(progressPanel.isCanceled())
				{
					progressPanel.appendLine("\tProcess interrupted. Deleting incomplete lines.");

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
                memoryAllocation = null;
            }
			progressPanel.appendLine("\tSEGY Microseismic Read");
			progressPanel.appendLine("\tMicroseismic Gather Written");
			progressPanel.appendLine("\tGather Header File Written");

			if(runTimer) totalTimer.stopPrint("Time to output all files:");
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackMicroseismic.readWritePreStackLines() failed.", e, StsException.WARNING);
			return false;
		}
	}	
}
