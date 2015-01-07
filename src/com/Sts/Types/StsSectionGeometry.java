
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsSectionGeometry extends StsSerialize
{
    boolean isVertical;
    byte geometryType = RIBBON;
    float rowF0, colF0, rowF1, colF1;
    int rowDirection = NONE; // if NONE, section is aligned on col or not at all
    int colDirection = NONE; // if NONE, section is aligned on row or not at all
    float dipAngle; // dip angle in degrees from vertical
    float dipDirection; // direction in degrees of dip (increasing depth)

    transient float newDipAngle;

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int ROWCOL = StsParameters.ROWCOL;

    static final int NONE = StsParameters.NONE;
    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    static final int RIGHT = StsParameters.RIGHT;
    static final int LEFT = StsParameters.LEFT;

    /** section consists of two lines with a twisted ribbon between */
    static public final byte RIBBON = 1;
    /** section consists of a series of more horizontal edges; each can be very different in shape, spacing, and number of points */
    static public final byte CURVED = 2;
    /** section consists of a series of ribs, typically fault-sticks. */
    static public final byte RIBS = 3;

    static final float nullValue = StsParameters.nullValue;

    static final float minGridDistance = 1.0f;

    static final long serialVersionUID = -4017646825327586445L;

    public StsSectionGeometry()
    {
    }
    public StsSectionGeometry(StsSection section, StsModel model)
    {
        this(section, getGeometryType(section), model);
    }
    public StsSectionGeometry(StsSection section, byte geometryType, StsModel model)
    {
        StsLine firstLine, lastLine;
        StsGridPoint firstGridPoint, lastGridPoint;

        this.geometryType = geometryType;

        firstLine = section.getFirstLine();
        lastLine = section.getLastLine();

        isVertical =  firstLine.getIsVertical() && lastLine.getIsVertical();

        StsGridDefinition gridDef = model.getGridDefinition();
        firstGridPoint = new StsGridPoint(firstLine.getTopPoint(), gridDef);
        lastGridPoint = new StsGridPoint(lastLine.getTopPoint(), gridDef);

        computeAlignment(firstGridPoint, lastGridPoint, false);
//        computeDipAngle(section);
    }

    public StsSectionGeometry(StsGridPoint gridPoint0, StsGridPoint gridPoint1, byte geometryType)
    {
        isVertical = true;
        this.geometryType = geometryType;
        computeAlignment(gridPoint0, gridPoint1, true);
    }

    static private byte getGeometryType(StsSection section)
    {
        int nSectionEdges = section.getNSectionEdges();
        if(nSectionEdges == 0)
            return RIBS;
        if(section.edgesAreStraight())
            return RIBBON;
        else
            return CURVED;
    }

    public void computeAlignment(StsGridPoint gridPoint0, StsGridPoint gridPoint1, boolean snap)
    {
        rowF0 = gridPoint0.rowF;
        colF0 = gridPoint0.colF;
        rowF1 = gridPoint1.rowF;
        colF1 = gridPoint1.colF;

        if(!isVertical || geometryType != RIBBON) return;

        float dRowF = rowF1 - rowF0;
        float dColF = colF1 - colF0;

        if(snap)
        {
            if(dRowF < minGridDistance && dRowF > -minGridDistance)
            {
                int row = StsMath.roundOffInteger(0.5f*(rowF0 + rowF1));
                rowF0 = (float)row;
                rowF1 = (float)row;
                gridPoint0.setRow(row);
                gridPoint1.setRow(row);
                dRowF = 0.0f;
            }
            if(dColF < minGridDistance && dColF > -minGridDistance)
            {
                int col = StsMath.roundOffInteger(0.5f*(colF0 + colF1));
                colF0 = (float)col;
                colF1 = (float)col;
                gridPoint0.setCol(col);
                gridPoint1.setCol(col);
                dColF = 0.0f;
            }
        }
        if(dRowF == 0.0f)
            colDirection = NONE;
        else if(dRowF > 0.0f)
            colDirection = PLUS;
        else
            colDirection = MINUS;

        if(dColF == 0.0f)
            rowDirection = NONE;
        else if(dColF > 0.0f)
            rowDirection = PLUS;
        else
            rowDirection = MINUS;
    }

    public void computeDipAngle(StsSection thisSection)
    {
        float[] dipDirectionVector = computeDipDirectionVector(RIGHT);
        dipDirection = StsMath.atan2(dipDirectionVector);

        if(isVertical)
            dipAngle = 0.0f;
        else
            dipAngle = thisSection.computeDipAngle(dipDirectionVector);
    }

    public float[] computeDipDirectionVector(int side)
    {
        float dipDirection;

        if(isRowColAligned())
        {
            int alignment = getRowOrCol();
            int direction = getRowColDirection();
            if(side == LEFT) direction = -direction;

            if(alignment == ROW)
            {
                if(direction == PLUS) return new float[] { 0.0f, -1.0f };
                else                  return new float[] { 0.0f, 1.0f };
            }
            else
            {
                if(direction == PLUS) return new float[] { 1.0f, 0.0f };
                else                  return new float[] { -1.0f, 0.0f };
            }
        }
        else
        {
            float dRowF, dColF;

            if(side == RIGHT)
            {
                dRowF = rowF1 - rowF0;
                dColF = colF1 - colF0;
            }
            else
            {
                dRowF = rowF0 - rowF1;
                dColF = colF0 - colF1;
            }

            float[] vector = new float[] { dRowF, -dColF };
            StsMath.normalizeVector(vector);
            return vector;
        }
    }

    public void dipAngleChanged(StsSection thisSection, float newDipAngle)
    {
        isVertical = false;
        this.newDipAngle = newDipAngle;
        thisSection.setInitialized(false);

        StsLine[] lines = thisSection.getLineList();
        int nLines = lines.length;
        for(int n = 0; n < nLines; n++)
        {
            lines[n].setIsVertical(false);
            lines[n].setInitialized(false);
        }
    }

    public void setNewDipAngle()
    {
        dipAngle = newDipAngle;
    }

    public boolean isCurved() { return geometryType == CURVED; }
    public boolean isVertical() { return isVertical; }
    public boolean isRibbon() { return geometryType == RIBBON; }
    public boolean isRibbed() { return geometryType == RIBS; }
    public boolean isPlanar() { return isVertical && geometryType == RIBBON; }
    public int getRowDirection() { return rowDirection; }
    public int getColDirection() { return colDirection; }
    public void setIsVertical(boolean isVertical) { this.isVertical = isVertical; }

    public float getDipAngle() { return dipAngle; }
    public float getNewDipAngle() { return newDipAngle; }
    public float getDipAngleChange() { return newDipAngle - dipAngle; }

    public boolean isRowColAligned()
    {
        int rowColAlignment = getRowOrCol();
        return rowColAlignment == ROW || rowColAlignment == COL;
    }

    // Section can be vertical or non-vertical.  If vertical return whether it is
    // aligned along ROW or COL.  If vertical, but not so aligned, return ROWCOL.
    // If nonVertical, return NONE.
    public int getRowOrCol()
    {
        if(geometryType != RIBBON || !isVertical) return NONE;

        if(rowDirection != NONE)
        {
            if(colDirection != NONE) return ROWCOL;
            else                  return ROW;
        }
        if(colDirection != NONE)     return COL;
        else                      return NONE;
    }

    public int getRowColDirection()
    {
        if(colDirection == 0.0f) return rowDirection;
        else if(rowDirection == 0.0f) return colDirection;
        else return NONE;
    }

    public byte getGeometryType()
    {
        if(colDirection == NONE)
        {
            if(rowDirection == PLUS)
                return StsSection.GEOM_ROW_PLUS;
            else if(rowDirection == MINUS)
                return StsSection.GEOM_ROW_MINUS;
        }
        else if(rowDirection == NONE)
        {
            if(colDirection == PLUS)
                return StsSection.GEOM_COL_PLUS;
            else if(colDirection == MINUS)
                return StsSection.GEOM_COL_MINUS;
        }
        return StsSection.GEOM_UNALIGNED;
    }
}
