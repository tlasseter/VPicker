package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Utilities.*;

import java.nio.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 11:16:33 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackFold extends StsDisplayPreStackAttribute
{
    boolean rangeSet = false;

    private StsPreStackFold(StsPreStackLineSet lineSet, StsModel model)
    {
        super(lineSet, model);
    }

    static public StsPreStackFold constructor(StsPreStackLineSet lineSet, StsModel model)
    {
        return new StsPreStackFold(lineSet, model);
    }

    public String getName() { return StsPreStackLineSet.ATTRIBUTE_FOLD_STRING; }

    public boolean computeBytes()
    {
        int maxFold = 0;
        int row, col, nTraces;
        int nRows = lineSet.nRows;
        int nCols = lineSet.nCols;
        StsPreStackLine[] lines = lineSet.lines;

        int[][] fold = new int[nRows][nCols];
        int nLines = lines.length;
        for (int n = 0; n < nLines; n++)
        {
            StsPreStackLine line = lines[n];
            if (line.isInline)
            {
                row = line.lineIndex;
                if (row < 0) row = nRows + line.lineIndex;
                int colMin = line.minGatherIndex;
                int colMax = line.maxGatherIndex;
                for (col = colMin; col <= colMax; col++)
                {
                    nTraces = line.getNTracesInGather(col);
                    maxFold = Math.max(maxFold, nTraces);
                    fold[row][col] = nTraces;
                }
            }
            else
            {
                col = line.lineIndex;
                int rowMin = line.minGatherIndex;
                int rowMax = line.maxGatherIndex;
                for (row = rowMin; row <= rowMax; row++)
                {
                    nTraces = line.getNTracesInGather(row);
                    maxFold = Math.max(maxFold, nTraces);
                    fold[row][col] = nTraces;
                }
            }
        }
        if (maxFold == 0) return false;

        if(colorscale == null)
            colorscale = getInitializeColorscale(StsSpectrumClass.SPECTRUM_RAINBOW, 0.0f, (float)maxFold);
        float scale = 254.0f/colorscale.getEditMax();

        byte[] byteFold = new byte[nRows * nCols];
        int n = 0;
        for (row = 0; row < nRows; row++)
        {
            for (col = 0; col < nCols; col++)
            {
                int traceFold = fold[row][col];
                if (traceFold == 0)
                {
                    byteFold[n++] = StsPreStackLineSet.nullByte;
                }
                else
                {
                    byteFold[n++] = StsMath.unsignedIntToUnsignedByte254( (int) (scale * fold[row][col] -1));
                }
            }
        }
        fold = null;
        byteBuffer = ByteBuffer.wrap(byteFold);
        return true;
    }
}
