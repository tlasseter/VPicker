
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Types.*;

public class StsSectionLinkGrid extends StsLinkGrid
{
    StsSection section;

    public StsSectionLinkGrid()
    {
    }

    public StsSectionLinkGrid(StsSection section, StsList edgeLinks)
    {
        super(section, edgeLinks);
        this.section = section;
    }
/*
    public void addInsideLinks()
    {
        for(int row = rowMin; row <= rowMax; row++)
            addInsideLinks(ROW, row);

        for(int col = colMin; col <= colMax; col++)
            addInsideLinks(COL, col);
    }

    public void addInsideLinks(int rowOrCol, int rowCol)
    {
        int linkIndex;
        StsGridLink prevLink, nextLink;
        int row, col;
        StsPoint point;
        StsGridSectionPoint gridPoint;
        StsGridLink newLink = null;

        try
        {
            nextLink = getFirstRowColLink(rowOrCol, rowCol);
            if(nextLink == null) return;

            linkIndex = StsGridLink.getLinkIndex(rowOrCol, PLUS);

            prevLink = nextLink;
            nextLink = prevLink.getLink(linkIndex);

            while(nextLink != null)
            {
                int start = StsMath.above(prevLink.getRowColF(rowOrCol));
                int end = StsMath.below(nextLink.getRowColF(rowOrCol));

                if(rowOrCol == ROW)
                {
                    row = rowCol;

                    for(col = start; col <= end; col++)
                    {
                        newLink = createNewLink(row, col);
                        newLink.connectToLinks(rowOrCol, prevLink, nextLink, grid);
                        addGridCrossingLinks(rowOrCol, prevLink, newLink);
                        prevLink = newLink;
                    }
                    if(newLink != null) addGridCrossingLinks(rowOrCol, newLink, nextLink);
                }
                else if(rowOrCol == COL)
                {
                    col = rowCol;

                    for(row = start; row <= end; row++)
                    {
                        newLink = createNewLink(row, col);
                        newLink.connectToLinks(rowOrCol, prevLink, nextLink, grid);
                        addGridCrossingLinks(rowOrCol, prevLink, newLink);
                        prevLink = newLink;
                    }
                    if(newLink != null) addGridCrossingLinks(rowOrCol, newLink, nextLink);
                }
                prevLink = nextLink;
                nextLink = prevLink.getLink(linkIndex);
            }
            StsException.systemError("StsLinkGrid.addInsidePoints() failed. " +
                "Could not create inside points on " + StsParameters.rowCol(rowOrCol) + rowCol);
        }
        catch(Exception e)
        {
            StsException.outputException("StsLinkGrid.addInsidePoints() failed. " +
                "Could not create inside points on " + StsParameters.rowCol(rowOrCol) + rowCol,
                e, StsException.WARNING);
        }
    }

    public StsGridLink createNewLink(int row, int col)
    {
        StsPoint point = grid.getPoint(row, col);
        StsGridSectionPoint gridPoint = new StsGridSectionPoint(point, (float)row, (float)col, null, section, false);
        return StsGridLink.constructInteriorLink(section, gridPoint);
    }

    private void addGridCrossingLinks(int rowOrCol, StsGridLink link0, StsGridLink link1)
    {
        StsGridLink prevLink, nextLink;
        StsGridSectionPoint gridPoint0, gridPoint1;

        gridPoint0 = (StsGridSectionPoint)link0.getPoint();
        gridPoint1 = (StsGridSectionPoint)link1.getPoint();

        StsList crossingPoints = StsGridSectionPoint.getGridCrossings(gridPoint0, gridPoint1, false);

        int nPnts = crossingPoints.getSize();
        if(nPnts <= 0) return;

        prevLink = link0;
        for(int n = 0; n < nPnts; n++)
        {
            StsGridSectionPoint gridPoint = (StsGridSectionPoint)crossingPoints.getElement(n);
            nextLink = new StsGridLink(section, gridPoint);
            prevLink.connectToNextLink(rowOrCol, nextLink, grid);
            prevLink = nextLink;
        }
        prevLink.connectToNextLink(rowOrCol, link1, grid);
    }
*/
}
