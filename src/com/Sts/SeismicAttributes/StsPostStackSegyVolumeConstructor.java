package com.Sts.SeismicAttributes;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version 1.0
 */

/** This is just a workhorse class as it's created and used just for computing a postStack attribute volume */
public class StsPostStackSegyVolumeConstructor extends StsVolumeConstructor
{
    StsSegyVolume segyVolume;
    StsSegyIO segyIO = null;
    StsCroppedBoundingBox cropBox;
    StsSegyData segyData;
    boolean readForward;

    int inputRowInc;
    int inputColInc;
    int inputSliceInc;

    int nInputFirstRow;
    int nInputLastRow;

    int nInputBytesPerBlock;

    int blockRowMin;
    int blockRowMax;

//    int nOutputBytesPerRow;
    int bytesPerTrace;
    boolean isCropped = false;
    /** number of attributes saved per traceHeader */
    int nAttributes;
    StsSegyVolume.SeismicLine[] lines;
    FileChannel inputSegyChannel;
    MappedByteBuffer inputBuffer = null;

    /** Indicates cube is regular (rectangular) */
    public boolean isRegular = true;

    /** attributes for each trace; sequential for each attribute [nAttributes][nTraces] */
    double[][] traceAttributes;
    /** number of traces written to attribute array */
    int nAttributeTrace = 0;
    /** Attribute records for extraction */
    StsSEGYFormatRec[] attributeRecords = null;

    static public byte[] byteTransparentTrace = null;
    static public float[] floatTransparentTrace = null;

    static byte[] segyTraceBytes;
    static byte[] traceBytes;
    static byte[] paddedTraceBytes;
    static float[] paddedTraceFloats;

    public boolean isOutputDataFloat() { return true; }

    public boolean initializeBlockInput() { return true; }

    public StsPostStackSegyVolumeConstructor(StsModel model, StsSegyVolume segyVolume, StsSeismicVolume volume, StsProgressDialog dialog)
    {
        initialize(model, segyVolume, dialog.getProgressPanel());
        this.dialog = dialog;
        createOutputVolume();
    }

    public StsPostStackSegyVolumeConstructor(StsModel model, StsSegyVolume segyVolume, StsProgressPanel panel)
    {
        initialize(model, segyVolume, panel);
        createOutputVolume();
    }

    private void initialize(StsModel model, StsSegyVolume segyVolume, StsProgressPanel panel)
    {
        this.model = model;
        cropBox = segyVolume.cropBox;
        outputVolume = new StsSeismicVolume(model, segyVolume, cropBox, "rw");
        this.panel = panel;
        attributeRecords = segyVolume.getSegyFormat().getRequiredTraceRecords();

        this.segyVolume = segyVolume;
        nInputRows = segyVolume.nRows;
        nInputCols = segyVolume.nCols;
        nInputSlices = segyVolume.nSlices;

        nInputFirstRow = cropBox.rowMin;
        nInputLastRow = cropBox.rowMax;
        isCropped = cropBox.isCropped;

        segyData = segyVolume.segyData;
        segyIO = segyVolume.segyIO;

        bytesPerTrace = segyData.bytesPerTrace;

        readForward = segyVolume.isXLineCCW;
        isRegular = segyVolume.isRegular;
        if (!isRegular)
        {
            lines = segyVolume.lines;
            paddedTraceBytes = new byte[cropBox.nSlices];
            paddedTraceFloats = new float[cropBox.nSlices];
            for (int n = 0; n < cropBox.nSlices; n++)
            {
                paddedTraceBytes[n] = -1;
                paddedTraceFloats[n] = StsParameters.nullValue;
            }
        }

        segyVolume.initializeAttributes();

        segyTraceBytes = new byte[bytesPerTrace];
        traceBytes = new byte[cropBox.nSlices];

        checkTransparentTrace(nInputSlices);

        attributeRecords = segyVolume.getSegyFormat().getRequiredTraceRecords();
        nAttributes = attributeRecords.length;
    }

    static public void checkTransparentTrace(int nSlices)
    {
        if (byteTransparentTrace != null && byteTransparentTrace.length > nSlices) return;
        byteTransparentTrace = new byte[nSlices];
        floatTransparentTrace = new float[nSlices];
        for (int n = 0; n < nSlices; n++)
        {
            byteTransparentTrace[n] = -1;
            floatTransparentTrace[n] = 0;
        }
        // TODO: classInitialize floatTransparentTrace to appropriate value
    }

    static public void computeVolumeWithDialog(StsModel model, StsSegyVolume segyVolume, StsSeismicVolume volume)
    {
        StsProgressDialog dialog = new StsProgressDialog(model.win3d, "Compute Seismic PostStack3d ", false);
        new StsPostStackSegyVolumeConstructor(model, segyVolume, volume, dialog);
    }

    protected String getFullStemname(StsModel model)
    {
        return model.getProject().getName() + "." + segyVolume.stemname;
    }

    public boolean allocateMemory()
    {
        memoryAllocation = StsMemAllocVolumeProcess.constructor(segyVolume, cropBox, usingOutputMappedBuffers);
        if (memoryAllocation == null) return false;
        nInputRowsPerBlock = memoryAllocation.nInputRowsPerBlock;
        nInputBytesPerBlock = memoryAllocation.nInputBytesPerBlock;
        nInputBlockTraces = memoryAllocation.nInputBlockTraces;
        nOutputRowsPerBlock = memoryAllocation.nOutputRowsPerBlock;
        nOutputSamplesPerInputBlock = memoryAllocation.nOutputSamplesPerInputBlock;
//        nOutputBytesPerRow = memoryAllocation.nOutputBytesPerRow;
        nOutputSamplesPerRow = memoryAllocation.nOutputSamplesPerRow;
        return true;
    }

    public boolean initializeBlockInput(int nBlock)
    {
        boolean isLastBlock = false;
        if (readForward)
        {
            nInputBlockFirstRow = nInputBlockLastRow + 1;
            nInputBlockLastRow += nInputRowsPerBlock;

            if (nInputBlockLastRow > nInputLastRow)
            {
                nInputBlockLastRow = nInputLastRow;
                isLastBlock = true;
            }
        }
        else
        {
            nInputBlockLastRow = nInputBlockFirstRow - 1;
            nInputBlockFirstRow -= nInputRowsPerBlock;

            if (nInputBlockFirstRow < nInputFirstRow) // for last block, possibly reduce block size
            {
                nInputBlockFirstRow = nInputFirstRow;
                inputPosition = segyData.fileHeaderSize + nInputFirstRow*nInputCols*bytesPerTrace;
                isLastBlock = true;
            }
        }
        if(isLastBlock)
        {
            nInputRowsPerBlock = nInputBlockLastRow - nInputBlockFirstRow + 1;
            nInputBlockTraces = nInputRowsPerBlock * nInputCols;
            nInputBytesPerBlock = nInputBlockTraces * bytesPerTrace;
            nOutputSamplesPerInputBlock = nInputRowsPerBlock*nOutputSamplesPerRow;
        }
        cropBox.adjustCropBlockRows(nInputBlockFirstRow, nInputBlockLastRow, nInputFirstRow);
        blockRowMin = cropBox.blockRowMin;
        blockRowMax = cropBox.blockRowMax;
        nOutputSamplesPerInputBlock = nOutputCols*nOutputSlices*(1 + (blockRowMax - blockRowMin)/inputRowInc);
        checkAttributesArray();
        return true;
    }

/*
    public void adjustCropBlockRows()
    {
        blockRowMin = StsMath.intervalRoundUp(nInputBlockFirstRow, inputRowInc);
        blockRowMax = StsMath.intervalRoundDown(nInputBlockLastRow, inputRowInc);

        nBlockRows = Math.max(0, (blockRowMax - blockRowMin) / inputRowInc + 1);
        nBlockSamples = nBlockRows * nSamplesPerRow;
        nBlockColSamples = nBlockRows * nCroppedSlices;
        nBlockSliceSamples = nBlockRows * nCols;

        outputRowMin = outputRowMax + 1;
        outputRowMax = outputRowMin + nBlockRows - 1;
    }
*/
    public  boolean processBlockInput(int nBlock, String mode ) { return processBlockInput(nBlock); }
    public boolean processBlockInput(int nBlock)
    {
 		if (runTimer) timer.getBlockInputTimer.start();

        try
        {
            if (isRegular)
            {
                if (debug)
                    System.out.println("Remapping segy file channel. inputPosition: " + inputPosition + " nInputBytes: " + nInputBytesPerBlock);
                try
                {
                    inputBuffer = inputSegyChannel.map(FileChannel.MapMode.READ_ONLY, inputPosition, nInputBytesPerBlock);
                }
                catch (IOException ioe)
                {
                    StsException.outputWarningException(this, "doProcessBlock", "inputSegyChannel.map() failed.", ioe);
                    return false;
                }
                if (readForward)
                {
                    inputPosition += nInputBytesPerBlock;
                }
                else
                {
                    inputPosition -= nInputBytesPerBlock;
                    if (inputPosition < segyData.fileHeaderSize)
                        inputPosition = segyData.fileHeaderSize;
                }

                if (debug) System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " nOutputFloatBytesPerInputBlock: " + nOutputSamplesPerInputBlock);
                outputFloatBuffer.map(outputPosition, nOutputSamplesPerInputBlock);
                outputPosition += nOutputSamplesPerInputBlock;

                if (runTimer)
                    StsVolumeConstructorTimer.getBlockInputTimer.stopAccumulateIncrementCountPrintInterval("   mapping of input and inline blocks for block " + nBlock + ":");

                if (!isCropped)
                {
                    if (!processBlock(nBlock))
                        return false;
                }
                else // isCropped
                {
                    if (!processCroppedBlock(outputFloatBuffer))
                        return false;
                }
            }
            else // not regular
            {
                if (!isCropped)
                {
                    if (!processIrregularBlock(nBlock))
                        return false;
                }
                else // isCropped
                {
                    if (!processCroppedIrregularBlock(nBlock))
                        return false;
                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "doProcessBlock", e);
            return false;
        }
        finally
        {
            clearMappedBuffers(nBlock);
        }
    }

    private void checkAttributesArray()
    {
        if (traceAttributes == null || traceAttributes[0].length != cropBox.nBlockSliceSamples)
            traceAttributes = new double[nAttributes][cropBox.nBlockSliceSamples];
    }

    private void clearMappedBuffers(int nBlock)
    {
        if(debug) StsException.systemDebug(this, "clearMappedBuffers", " input and outputFloat buffers for block: " + nBlock);

        if (inputBuffer != null)
        {
            if (runTimer) StsVolumeConstructorTimer.getClearInputBlockTimer.start();
            inputBuffer.clear();
            StsToolkit.clean(inputBuffer);
            if (runTimer) StsVolumeConstructorTimer.getClearInputBlockTimer.stopPrint("Time to clear input mapped buffers for block " + nBlock + ":");
        }

        if(outputFloatBuffer != null)
        {
//            outputFloatBuffer.force();
            if(runTimer)
                outputFloatBuffer.clearDebug("input mapped buffers for block " + nBlock, StsVolumeConstructorTimer.getClearInputBlockTimer);
            else
                outputFloatBuffer.clear();
        }
    }

    private boolean processBlock(int nBlock)
    {
        int nBlockTrace = -1;
        int bytesPerTrace = segyData.bytesPerTrace;

        try
        {
            if (runTimer) StsVolumeConstructorTimer.getBlockInputTimer.start();
            if (readForward)
            {
                inputBuffer.position(0);
                for (nBlockTrace = 0; nBlockTrace < nInputBlockTraces; nBlockTrace++)
                {
                    inputBuffer.get(segyTraceBytes);
                    extractHeaderAttributes(segyTraceBytes, nBlockTrace);
                    segyIO.processTrace(segyTraceBytes, outputFloatBuffer, nBlockTrace);
                }
            }
            else // backwards
            {
                for (int row = cropBox.nBlockRows - 1; row >= 0; row--)
                {
                    nBlockTrace = row * nInputCols;
                    for (int col = 0; col < nInputCols; col++, nBlockTrace++)
                    {
                        int pos = nBlockTrace * bytesPerTrace;
                        inputBuffer.position(pos);
                        inputBuffer.get(segyTraceBytes);
                        extractHeaderAttributes(segyTraceBytes, nBlockTrace);
                        segyIO.processTrace(segyTraceBytes, outputFloatBuffer, nBlockTrace);
                    }
                }
            }
            if (runTimer)
                StsVolumeConstructorTimer.getBlockInputTimer.stopAccumulateIncrementCountPrintInterval("process " + nInputBlockTraces + " traces for block " + nBlock + ":");
                return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "processBlock", "failed at  trace " + nBlockTrace, e);
            return false;
        }
    }

    private boolean processCroppedBlock(StsMappedFloatBuffer outputRowFloatBuffer)
    {
        int row = 0, col = 0, n = 0, pos = 0;

        try
        {
            int rowMin = blockRowMin - nInputBlockFirstRow;
            int rowMax = blockRowMax - nInputBlockFirstRow;
            int rowInc = inputRowInc;
            int colMin = cropBox.colMin;
            int colMax = cropBox.colMax;
            int colInc = cropBox.colInc;
            int sliceMin = cropBox.sliceMin;
            int sliceMax = cropBox.sliceMax;
            int sliceInc = cropBox.sliceInc;

            int nOutputTrace = 0;
            if(readForward)
            {
                for(row = rowMin; row <= rowMax; row+=rowInc, nOutputTrace++)
                {
                    for (col = colMin; col <= colMax; col += colInc)
                    {
                        int nTrace = row * nInputCols + col;
                        pos = nTrace * segyData.bytesPerTrace;
                        inputBuffer.position(pos);
                        inputBuffer.get(segyTraceBytes);
                        extractHeaderAttributes(segyTraceBytes, nOutputTrace);
                        processCroppedTrace(segyTraceBytes, outputRowFloatBuffer, sliceMin, sliceMax, sliceInc, nTrace);
                    }
                }
            }
            else // readBackwards
            {
                for(row = rowMax; row >= rowMin; row-=rowInc, nOutputTrace++)
                {
                    for (col = colMin; col <= colMax; col += colInc)
                    {
                        int nTrace = row * nInputCols + col;
                        pos = nTrace * segyData.bytesPerTrace;
                        inputBuffer.position(pos);
                        inputBuffer.get(segyTraceBytes);
                        extractHeaderAttributes(segyTraceBytes, nOutputTrace);
                        processCroppedTrace(segyTraceBytes, outputRowFloatBuffer, sliceMin, sliceMax, sliceInc, nTrace);
                    }
                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "processCroppedBlock",
                    "failed at row " + row + " col " + col + " sample " + n + " pos " + pos, e);
            return false;
        }
    }

    final private boolean processCroppedTrace(byte[] segyTraceBytes, StsMappedFloatBuffer rowFloatBuffer, int sliceMin, int sliceMax, int sliceInc, int nTrace)
    {
        int n = 0, pos = 0, s = 0;
        int bytesPerSample = segyData.bytesPerSample;
        try
        {
            int offset = segyData.traceHeaderSize + sliceMin * bytesPerSample;
            //Process Samples
            for (s = sliceMin; s <= sliceMax; s += sliceInc)
            {
                segyIO.processSampleBytes(segyTraceBytes, offset, rowFloatBuffer);
                offset += sliceInc * bytesPerSample;
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "processCroppedTrace",
                    " failed at  trace " + nTrace + " sample " + s + " pos " + pos, e);
            return false;
        }
    }

    private boolean processBlockX(int nBlock)
    {
        int nBlockTrace = -1;
        int bytesPerTrace = segyData.bytesPerTrace;

        try
        {
            if (runTimer) StsVolumeConstructorTimer.getBlockInputTimer.start();
            if (readForward)
            {
                inputBuffer.position(0);
                for (nBlockTrace = 0; nBlockTrace < nInputBlockTraces; nBlockTrace++)
                {
                    inputBuffer.get(segyTraceBytes);
                    extractHeaderAttributes(segyTraceBytes, nBlockTrace);
                    segyIO.processTrace(segyTraceBytes, outputFloatBuffer, nBlockTrace);
                }
            }
            else // backwards
            {
                for (int row = cropBox.nBlockRows - 1; row >= 0; row--)
                {
                    nBlockTrace = row * nInputCols;
                    for (int col = 0; col < nInputCols; col++, nBlockTrace++)
                    {
                        int pos = nBlockTrace * bytesPerTrace;
                        inputBuffer.position(pos);
                        inputBuffer.get(segyTraceBytes);
                        extractHeaderAttributes(segyTraceBytes, nBlockTrace);
                        segyIO.processTrace(segyTraceBytes, outputFloatBuffer, nBlockTrace);
                    }
                }
            }
            if (runTimer)
                StsVolumeConstructorTimer.getBlockInputTimer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");

            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "processBlock", "failed at  trace " + nBlockTrace, e);
            return false;
        }
    }

    private boolean processIrregularBlock(int nBlock)
    {
        StsSEGYFormat.TraceHeader firstBlockTraceHeader = null, lastBlockTraceHeader = null;
        StsSEGYFormat.TraceHeader firstRowTraceHeader = null, lastRowTraceHeader = null;
        int bytesPerTrace = segyData.bytesPerTrace;
        boolean isLittleEndian = segyData.isLittleEndian;
        byte[] prevTraceBytes = new byte[bytesPerTrace]; // prev segy trace bytes
        byte[] nextTraceBytes = new byte[bytesPerTrace]; // next segy trace bytes
        double[] nextTraceAttributes = new double[nAttributes];
        double[] prevTraceAttributes = new double[nAttributes];
        int attributeIndexXline = segyVolume.getAttributeIndex(StsSEGYFormat.XLINE_NO);
        int attributeIndexNTrace = segyVolume.getAttributeIndex(StsSEGYFormat.TRACENO);
        int row = -1;
        int col = -1;
        int n = -1;
        int nOutputBytes = 0;
        int nRowTracesWritten = 0;
        try
        {
//            if(panel != null)
//                panel.setCurrentActionLabel("Processing irregular SegY file block " + nBlock);
            if (runTimer) StsVolumeConstructorTimer.getBlockInputTimer.start();

            row = nInputBlockFirstRow;
            while (firstBlockTraceHeader == null && row <= nInputBlockLastRow)
                firstBlockTraceHeader = lines[row++].firstTraceHeader;
            if (firstBlockTraceHeader == null) return false;
            int nFirstBlockTrace = firstBlockTraceHeader.nTrace;

            row = nInputBlockLastRow;
            while (lastBlockTraceHeader == null && row >= nInputBlockFirstRow)
                lastBlockTraceHeader = lines[row--].lastTraceHeader;
            if (lastBlockTraceHeader == null) return false;
            int nLastBlockTrace = lastBlockTraceHeader.nTrace;

            long blockInputPosition = (long) segyData.fileHeaderSize + (long) nFirstBlockTrace * (long) bytesPerTrace;
            int nBlockTraces = (nLastBlockTrace - nFirstBlockTrace + 1);
            nInputBytesPerBlock = nBlockTraces * segyData.bytesPerTrace;

            if(debug) System.out.println("Remapping segy file channel. blockInputPosition: " + blockInputPosition + " nInputBytes: " + nInputBytesPerBlock);
            inputBuffer = inputSegyChannel.map(FileChannel.MapMode.READ_ONLY, blockInputPosition, nInputBytesPerBlock);

            if(debug) System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " cropBox.nBlockSamples: " + cropBox.nBlockSamples);
            outputFloatBuffer.map(outputPosition, cropBox.nBlockSamples);
            outputPosition += cropBox.nBlockSamples;

            if (runTimer)
            {
                StsVolumeConstructorTimer.getBlockInputTimer.stopPrint("NIO mapping of input and inline blocks for block " + nBlock + ":");
                StsVolumeConstructorTimer.getBlockInputTimer.start();
            }

            // loop for indexed lines

            int nFirstRow, nLastRow, rowInc;
            if(readForward)
            {
                nFirstRow = nInputBlockFirstRow;
                nLastRow = nInputBlockLastRow;
                rowInc = 1;
            }
            else
            {
                nFirstRow = nInputBlockLastRow;
                nLastRow = nInputBlockFirstRow;
                rowInc = -1;
            }
            nOutputBytes = 0;
            row = nFirstRow;
            while (true)
            {
                nRowTracesWritten = 0;
                col = 0;

                firstRowTraceHeader = lines[row].firstTraceHeader;
                lastRowTraceHeader = lines[row].lastTraceHeader;

                if (firstRowTraceHeader == null && lastRowTraceHeader == null)
                {
                    int nPaddedStartTraces = nInputCols;

                    for (n = 0; n < nPaddedStartTraces; n++)
                    {
                        col++;
                        outputFloatBuffer.put(paddedTraceFloats);
                        nOutputBytes += nInputSlices;
                    }
                    nRowTracesWritten += nPaddedStartTraces;
                    if (row == nLastRow) break;
                    row += rowInc;
                    continue;
                }

                int nFirstRowTrace = firstRowTraceHeader.nTrace;
                int nLastRowTrace = lastRowTraceHeader.nTrace;
                int nBlockTrace;
                if (readForward)
                    nBlockTrace = (row - nInputBlockFirstRow) * nInputCols;
                else
                    nBlockTrace = (nInputBlockLastRow - row) * nInputCols;

                // write missing traces at start of line
                int nPaddedStartTraces = getTraceColIndex(firstRowTraceHeader);
                if(nPaddedStartTraces > 0)
                {
                    panel.appendLine("Padding initial traces for inline: " + firstRowTraceHeader.iLine +
                                " crosslines: " + segyVolume.rowNumInc + "-" + firstRowTraceHeader.xLine);
                }
                for (n = 0; n < nPaddedStartTraces; n++)
                {
                    col++;
                    outputFloatBuffer.put(paddedTraceFloats);
                    nOutputBytes += nInputSlices;
                }
                nRowTracesWritten += nPaddedStartTraces;
                nBlockTrace += nPaddedStartTraces;
                // write traces and padded traces in between
                int pos = (nFirstRowTrace - nFirstBlockTrace) * bytesPerTrace;
                inputBuffer.position(pos);
                inputBuffer.get(nextTraceBytes);
                // we are assuming first and last trace on like are ok header-attribute-wise
                // if an intermediate trace has a bad header, we will simply pad it out
                // if any intermediate traces are missing, we will pad those out as well
                extractHeaderAttributes(nextTraceBytes, nextTraceAttributes);
                boolean nextTraceOk = headerAttributesOk(nextTraceAttributes, attributeIndexXline, firstRowTraceHeader, lastRowTraceHeader);
                float nextXLine = (float)nextTraceAttributes[attributeIndexXline];
                float prevXLine;
                boolean prevTraceOk;
                float lastXLine = segyVolume.colNumMin - segyVolume.colNumInc;
                for (n = nFirstRowTrace; n < nLastRowTrace; n++)
                {
                    System.arraycopy(nextTraceBytes, 0, prevTraceBytes, 0, bytesPerTrace);
                    System.arraycopy(nextTraceAttributes, 0, prevTraceAttributes, 0, nAttributes);
                    prevXLine = nextXLine;
                    prevTraceOk = nextTraceOk;
                    int nMissingTraces = 0;
                    // loop until we find a nextTrace whose header is ok and nextXLine is sequenced properly from prevXLine (nMissingTraces >= 0)
                    while(n < nLastRowTrace)
                    {
                        inputBuffer.get(nextTraceBytes);
                        extractHeaderAttributes(nextTraceBytes, nextTraceAttributes);
                        // n++;
                        nextTraceOk = headerAttributesOk(nextTraceAttributes, attributeIndexXline, firstRowTraceHeader, lastRowTraceHeader);
                        if(nextTraceOk)
                        {
                            nextXLine = (float)nextTraceAttributes[attributeIndexXline];
                            // If nMissing traces == 0, nextXLine is in proper sequence with no missing trace between prevTrace and nextTrace.
                            // If nMissing traces > 0, we have nMissingTraces between prevTrace and nextTrace.
                            // If nMissing traces < 0, we have a xLine numbering problem: nextXLine == prevXLine or has reversed direction in xLine sequence.
                            //    in this case, we will subsequently output the prevTrace, but not the nextTrace.
                            nMissingTraces = Math.round((nextXLine - prevXLine)/segyVolume.colNumInc) - 1;
                            if(nMissingTraces >= 0) break;
                            // if nMissingTraces < 0, we have repeat trace XLine number or reverse-direction XLine number
                            nextTraceOk = false;
                            int nTrace = (int)nextTraceAttributes[attributeIndexNTrace];
                            int nPrevTrace = nTrace - 1;
                            panel.appendLine("XLine numbering wrong for traces: " + nPrevTrace + " [" + prevXLine + "] and " + nTrace + " [" + nextXLine + "].");
                            n++;
                        }
                    }

                    // write prevTrace
                    if(prevTraceOk)
                    {
                        putHeaderAttributes(nBlockTrace, prevTraceAttributes);
                        if (!segyIO.processTrace(prevTraceBytes, outputFloatBuffer, n))
                        {
                            StsException.systemError("Failed to process trace: " + n + " at row: " + row + " col: " + col + "\n" +
                                    "    nFirstBlockTrace: " + nFirstBlockTrace + " nLastBlockTrace: " + nLastBlockTrace + " nOutputBytes: " + nOutputBytes);
                            return false;
                        }
                        col++;
                        nOutputBytes += nInputSlices;
                        nRowTracesWritten++;
                        nBlockTrace++;
                        lastXLine = prevXLine;
                    }
                    // between prevTrace and nextTrace, write in missing traces
                    // nextTraceOk and nextTraceOk must be true if nMissingTraces > 0
                    if(nMissingTraces > 0)
                    {
                        float nPaddedFirstXLine = prevXLine + 1;
                        float nPaddedLastXLine = nextXLine - 1;
                        if(nextTraceOk)
                        {
                            panel.appendLine("Padding traces for inline: " + firstRowTraceHeader.iLine +
                                    " crosslines: " + nPaddedFirstXLine + "-" + nPaddedLastXLine);
                        }
                        else
                        {
                            panel.appendLine("Padding bad traces for inline: " + firstRowTraceHeader.iLine +
                                    " crosslines: " + nPaddedFirstXLine + "-" + nPaddedLastXLine);
                        }

                        for (int i = 0; i < nMissingTraces; i++)
                        {
                            col++;
    //							inlineBuffer.put(paddedTraceBytes);
                            if (outputFloatBuffer != null)
                                outputFloatBuffer.put(paddedTraceFloats);
                            nOutputBytes += nInputSlices;
                        }
                        nRowTracesWritten += nMissingTraces;
                        nBlockTrace += nMissingTraces;
                    }
                }
                // write last trace
                if(nextTraceOk)
                {
                    putHeaderAttributes(nBlockTrace, nextTraceAttributes);
                    if (!segyIO.processTrace(nextTraceBytes, outputFloatBuffer, nLastRowTrace))
                    {
                        StsException.systemError("Failed to process trace: " + n + " at row: " + row + " col: " + col);
                        return false;
                    }
                    col++;
                    nOutputBytes += nInputSlices;
                    nRowTracesWritten++;
                    nBlockTrace++;
                    lastXLine = nextXLine;
                }
                // write missing traces at end of line
                float nPaddedFirstXLine = lastXLine + segyVolume.colNumInc;
                float nPaddedLastXLine = segyVolume.colNumMax;
                int nPaddedEndTraces = Math.round((nPaddedLastXLine - nPaddedFirstXLine) / segyVolume.colNumInc) + 1;
                if(nPaddedEndTraces > 0)
                {
                    if(nPaddedEndTraces == 1)
                        panel.appendLine("Padding trace for inline: " + firstRowTraceHeader.iLine + " crossline: " + nPaddedFirstXLine);
                    else
                        panel.appendLine("Padding traces for inline: " + firstRowTraceHeader.iLine +
                                " crosslines: " + nPaddedFirstXLine + "-" + nPaddedLastXLine);
                }

                for (int i = 0; i < nPaddedEndTraces; i++)
                {
                    col++;
                    outputFloatBuffer.put(paddedTraceFloats);
                    nOutputBytes += nInputSlices;
                }
                nRowTracesWritten += nPaddedEndTraces;
                nBlockTrace += nPaddedEndTraces;

                if (nRowTracesWritten != nInputCols)
                {
                    StsException.systemError("readIrregularBlock() error. Row: " + row +
                            " Wrote " + nRowTracesWritten + ". Should be: " + nInputCols + ".\n" +
                            " first trace: " + firstRowTraceHeader.toString() + " last trace: " +
                            lastRowTraceHeader.toString() + ".\n");
                    return false;
                }
                if (row == nLastRow)
                    break;
                row += rowInc;
            }
            if (nOutputBytes != cropBox.nBlockSamples)
            {
                StsException.systemError("processIrregularBlock() error. nOutputBytes: " + nOutputBytes + ". Should be: " +
                        cropBox.nBlockSamples);
                return false;
            }

            if (runTimer)
                StsVolumeConstructorTimer.getBlockInputTimer.stopPrint("Block " + nBlock + " processed " + nBlockTraces +
                                                      " irregular block traces into floats and write of inline block:");

            return true;
        }
        catch (Exception e)
        {
            int nSamplesPerRow = nInputCols * nInputSlices;
            int nRowStartBytes = (row - nInputBlockFirstRow) * nSamplesPerRow;
            int nRowEndBytes = nRowStartBytes + nSamplesPerRow;
            StsException.outputException("StsSegyVolume.readIrregularSegYFile failed while writing line index: " +
                    row + " xLine index: " + col, e, StsException.WARNING);
            String message = new String("Failed to load volume " + segyVolume.stemname +
                    " at line index: " + row + " at xLine index: " + col + ".\n" +
                    " nBytes written: " + nOutputBytes + ". Should be between " + nRowStartBytes +
                    " and " + nRowEndBytes + ".\n" +
                    " first trace: " + firstRowTraceHeader.toString() + " last trace: " +
                    lastRowTraceHeader.toString() + ".\n" +
                    "Error: " + e.getMessage());
            new StsMessage(model.win3d, StsMessage.WARNING, message);
            return false;
        }
    }

    private String traceDescription(double[] traceAttributes)
    {
        int indexXline = segyVolume.getAttributeIndex(StsSEGYFormat.XLINE_NO);
        int indexIline = segyVolume.getAttributeIndex(StsSEGYFormat.ILINE_NO);
        int indexNTrace = segyVolume.getAttributeIndex(StsSEGYFormat.XLINE_NO);
        return "trace " + traceAttributes[indexNTrace] + " iline " + traceAttributes[indexIline]  + " xline " + traceAttributes[indexXline];
    }

    private boolean headerAttributesOk(double[] nextTraceAttributes, int attributeIndexXline,
                                       StsSEGYFormat.TraceHeader firstRowTraceHeader, StsSEGYFormat.TraceHeader lastRowTraceHeader)
    {
        float nextXLine = (float)nextTraceAttributes[attributeIndexXline];
        return StsMath.betweenInclusive(nextXLine, firstRowTraceHeader.xLine, lastRowTraceHeader.xLine);
    }

    private boolean processCroppedIrregularBlock(int nBlock)
    {
        if(cropBox.nBlockRows <= 0) return true;
        StsSEGYFormat.TraceHeader firstBlockTraceHeader = null, lastBlockTraceHeader = null;
        StsSEGYFormat.TraceHeader firstRowTraceHeader = null, lastRowTraceHeader = null;
        int bytesPerTrace = segyData.bytesPerTrace;
        int bytesPerSample = segyData.bytesPerSample;
        byte[] prevTraceBytes = new byte[bytesPerTrace]; // prev segy trace bytes
        byte[] nextTraceBytes = new byte[bytesPerTrace]; // next segy trace bytes

//        StsSEGYFormat.TraceHeader prevRowTraceHeader, nextRowTraceHeader;
        int row = -1;
        int col = -1;
        int n = -1;
        int pos = -1;
        int nOutputBytes = 0;
        int nRowTracesWritten = 0;

        int rowMin = cropBox.blockRowMin;
        int rowMax = cropBox.blockRowMax;
        int rowInc = inputRowInc;
        int colMin = cropBox.colMin;
        int colMax = cropBox.colMax;
        int colInc = cropBox.colInc;
        int sliceMin = cropBox.sliceMin;
        int sliceMax = cropBox.sliceMax;
        int sliceInc = cropBox.sliceInc;

        try
        {
            if (runTimer) StsVolumeConstructorTimer.getBlockInputTimer.start();
            //if (panel != null)
            //panel.setCurrentActionLabel("Processing irregular SegY file block " + nBlock);
            panel.appendLine("Processing irregular SegY file block " + nBlock);

            row = rowMin;
            while (firstBlockTraceHeader == null)
                firstBlockTraceHeader = lines[row++].firstTraceHeader;
            if (firstBlockTraceHeader == null)
                return false;
            int nFirstBlockTrace = firstBlockTraceHeader.nTrace;

            row = rowMax;
            while (lastBlockTraceHeader == null)
                lastBlockTraceHeader = lines[row--].lastTraceHeader;
            if (lastBlockTraceHeader == null)
                return false;
            int nLastBlockTrace = lastBlockTraceHeader.nTrace;

            inputPosition = (long) segyData.fileHeaderSize + (long) nFirstBlockTrace * (long) bytesPerTrace;
            int nBlockTraces = (nLastBlockTrace - nFirstBlockTrace + 1);
            nInputBytesPerBlock = nBlockTraces * bytesPerTrace;

            if (runTimer) StsVolumeConstructorTimer.getBlockInputTimer.start();

            if (debug) System.out.println("Remapping segy file channel. inputPosition: " + inputPosition + " nInputBytes: " + nInputBytesPerBlock);
            inputBuffer = inputSegyChannel.map(FileChannel.MapMode.READ_ONLY, inputPosition, nInputBytesPerBlock);

            if (debug) System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " cropBox.nBlockSamples: " + cropBox.nBlockSamples);
            outputFloatBuffer.map(outputPosition, cropBox.nBlockSamples);
            outputPosition += cropBox.nBlockSamples;
            if (runTimer)
            {
                StsVolumeConstructorTimer.getBlockInputTimer.stopPrint("NIO mapping of input and inline blocks for block " + nBlock + ":");
                StsVolumeConstructorTimer.getBlockInputTimer.start();
            }

            int nFirstRow, nLastRow;
            if(readForward)
            {
                nFirstRow = rowMin;
                nLastRow = rowMax;
            }
            else
            {
                nFirstRow = rowMax;
                nLastRow = rowMin;
                rowInc = -rowInc;
            }
            nOutputBytes = 0;
            row = nFirstRow;
            while (true)
            {
                nRowTracesWritten = 0;
                col = -1;

                firstRowTraceHeader = lines[row].firstTraceHeader;
                lastRowTraceHeader = lines[row].lastTraceHeader;

                if (firstRowTraceHeader == null && lastRowTraceHeader == null)
                {
                    int nPaddedStartTraces = cropBox.nCols;

                    for (n = 0; n < nPaddedStartTraces; n++)
                    {
                        col++;
//						inlineBuffer.put(paddedTraceBytes);
                        outputFloatBuffer.put(paddedTraceFloats);
                        nOutputBytes += cropBox.nSlices;
                    }
                    nRowTracesWritten += nPaddedStartTraces;
                    continue;
                }

                int nFirstTrace = firstRowTraceHeader.nTrace;
                int nLastTrace = lastRowTraceHeader.nTrace;

                int nNextTrace = nFirstTrace;
                pos = (nNextTrace - nFirstBlockTrace) * bytesPerTrace;
                inputBuffer.position(pos);
                inputBuffer.get(nextTraceBytes);
                int nextTraceCol = getTraceColIndex(nextTraceBytes);
                for (col = colMin; col <= colMax; col += colInc)
                {
                    // write missing trace
                    if (col < nextTraceCol)
                    {
                        outputFloatBuffer.put(paddedTraceFloats);
                        nOutputBytes += cropBox.nSlices;
                        nRowTracesWritten++;
                    }
                    // write trace
                    else if (col == nextTraceCol)
                    {
                        n = nNextTrace - nFirstBlockTrace;
                        int nOutputTrace = (int) ((row * (colMax - colMin)) + (col - colMin));
// Does not work yet                        extractHeaderAttributes(segyTraceBytes, nOutputTrace);
                        if (!processCroppedTrace(nextTraceBytes, outputFloatBuffer, sliceMin, sliceMax, sliceInc, n))
                        {
                            StsException.systemError("Failed to process trace: " + n + " at row: " + row + " col: " + col +
                                    "\n" +
                                    "    nFirstBlockTrace: " + nFirstBlockTrace + " nLastBlockTrace: " +
                                    nLastBlockTrace + " nOutputBytes: " + nOutputBytes);
                            return false;
                        }
                        nOutputBytes += cropBox.nSlices;
                        nRowTracesWritten++;
                    }
                    else // col > nextTraceCol
                    {
                        while (nextTraceCol != -1 && nextTraceCol < col)
                        {
                            nNextTrace++;
                            if (nNextTrace <= nLastTrace)
                            {
                                pos = (nNextTrace - nFirstBlockTrace) * bytesPerTrace;
                                inputBuffer.position(pos);
                                inputBuffer.get(nextTraceBytes);
                                nextTraceCol = getTraceColIndex(nextTraceBytes);
                            }
                            else // nNextTrace > nLastTrace: we're out of traces
                            {
                                nextTraceCol = -1;
                            }
                        }
                        if (nextTraceCol == col)
                        {
                            n = nNextTrace - nFirstBlockTrace;
                            int nOutputTrace = (int) ((row * (colMax - colMin)) + (col - colMin));
// Does not work yet                            extractHeaderAttributes(segyTraceBytes, nOutputTrace);
                            if (!processCroppedTrace(nextTraceBytes, outputFloatBuffer, sliceMin, sliceMax, sliceInc, n))
                            {
                                StsException.systemError("Failed to process trace: " + n + " at row: " + row + " col: " + col +
                                        "\n" +
                                        "    nFirstBlockTrace: " + nFirstBlockTrace + " nLastBlockTrace: " +
                                        nLastBlockTrace + " nOutputBytes: " + nOutputBytes);
                                return false;
                            }
                        }
                        else
                        {
                            outputFloatBuffer.put(paddedTraceFloats);
                        }

                        nOutputBytes += cropBox.nSlices;
                        nRowTracesWritten++;

                    }
                }
                if (nRowTracesWritten != cropBox.nCols)
                {
                    StsException.systemError("readIrregularBlock() error. Row: " + row +
                            " Wrote " + nRowTracesWritten + ". Should be: " + cropBox.nCols + ".\n" +
                            " first trace: " + firstRowTraceHeader.toString() + " last trace: " +
                            lastRowTraceHeader.toString() + ".\n");
                    return false;
                }
                if (row == nLastRow) break;
                row += rowInc;
            }
            if (nOutputBytes != cropBox.nBlockSamples)
            {
                StsException.systemError("processIrregularBlock() error. nOutputBytes: " + nOutputBytes + ". Should be: " + cropBox.nBlockSamples);
                return false;
            }

            if (runTimer)
                StsVolumeConstructorTimer.getBlockInputTimer.stopPrint("Block " + nBlock + " processed " + nBlockTraces +
                                                      " irregular cropped block traces into floats and write of inline block:");

            return true;
        }
        catch (Exception e)
        {
            int nSamplesPerRow = nInputCols * nInputSlices;
            int nRowStartBytes = (row - nInputBlockFirstRow) * nSamplesPerRow;
            int nRowEndBytes = nRowStartBytes + nSamplesPerRow;
            StsException.outputException("StsSegyVolume.readIrregularSegYFile failed while writing line index: " +
                    row + " xLine index: " + col, e, StsException.WARNING);
            String message = new String("Failed to load volume " + segyVolume.stemname +
                    " at line index: " + row + " at xLine index: " + col + ".\n" +
                    " nBytes written: " + nOutputBytes + ". Should be between " + nRowStartBytes +
                    " and " + nRowEndBytes + ".\n" +
                    " first trace: " + StsSEGYFormat.toStringStatic(firstRowTraceHeader) + " last trace: " +
                    StsSEGYFormat.toStringStatic(lastRowTraceHeader) + ".\n" +
                    "Error: " + e.getMessage());
            new StsMessage(model.win3d, StsMessage.WARNING, message);
            return false;
        }
    }

    private void extractHeaderAttributes(byte[] segyTraceBytes, int traceNum)
    {
        segyData.extractHeaderAttributes(segyTraceBytes, traceNum, attributeRecords, traceAttributes);
    }

    private void extractHeaderAttributes(byte[] segyTraceBytes, double[] traceAttributes)
    {
        segyData.extractHeaderAttributes(segyTraceBytes, attributeRecords, traceAttributes);
    }

    private void putHeaderAttributes(int traceNum, double[] headerAttributes)
    {
		for(int i = 0; i < nAttributes; i++)
			traceAttributes[i][traceNum] = headerAttributes[i];
    }

    private int getTraceColIndex(StsSEGYFormat.TraceHeader traceHeader)
    {
        return segyVolume.getTraceColIndex(traceHeader);
    }

    private int getTraceColIndex(byte[] segyTraceBytes)
    {
        return segyVolume.getTraceColIndex(segyTraceBytes);
    }

    public void initializeVolumeInput()
    {
        super.initializeVolumeInput();
        inputSegyChannel = segyData.randomAccessSegyFile.getChannel();
        outputFloatBuffer = outputVolume.createMappedFloatRowBuffer();

        if (readForward)
        {
            inputPosition = (long) segyData.fileHeaderSize + (long) cropBox.rowMin * (long) nInputCols * (long) bytesPerTrace;
            nInputBlockLastRow = cropBox.rowMin - 1;
        }
        else
        {
            inputPosition = (long) segyData.fileHeaderSize + (long) (cropBox.rowMax + 1) * (long) nInputCols * (long) bytesPerTrace - (long) nInputBytesPerBlock;
            nInputBlockFirstRow = cropBox.rowMax + 1;
        }
        outputPosition = 0;
        segyIO.initializeVolumeHistogramSamples();
    }

    protected void initializeGridRange()
    {
        nInputRows = memoryAllocation.nInputRows;
        nInputCols = memoryAllocation.nInputCols;
        nInputSlices = memoryAllocation.nInputSlices;
        nOutputRows = memoryAllocation.nOutputRows;
        nOutputCols = memoryAllocation.nOutputCols;
        nOutputSlices = memoryAllocation.nOutputSlices;
        inputRowInc = cropBox.rowInc;
        inputColInc = cropBox.colInc;
        inputSliceInc = cropBox.sliceInc;
    }

    public void initializeVolumeOutput()
    {
        outputVolume.setDataMin(segyIO.dataMin);
        outputVolume.setDataMax(segyIO.dataMax);
        super.initializeVolumeOutput();
    }

    public void finalizeVolumeOutput()
    {
//        outputVolume.checkFlipRowNumOrder();
        super.finalizeVolumeOutput();
    }
}
