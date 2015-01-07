
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Types.*;

import java.util.*;

public class QuadStrip
{
    public int rowNumber = -1;
    public int firstCol = -1;
    public int lastCol = -1;
    public int ID = -1;
    public byte cellType;
    public Object[][] polygons;

    public QuadStrip() { }

    public boolean fieldsAreValid()
    {
        if (rowNumber==-1) return false;
        if (firstCol==-1) return false;
        if (lastCol==-1) return false;
        return true;
    }

    public void addPolygons(StsList cellPolygons, int col)
    {
        if(cellPolygons == null) return;

        if(polygons == null)
        {
            int nCells = lastCol - firstCol;
            if(nCells <= 0) return;
            polygons = new Object[nCells][];
        }
        if(col < firstCol || col > lastCol)
        {
            StsException.systemError("QuadStrip.addPolygon() failed." +
                " col: " + col + " is out of range " + firstCol + " -  " + lastCol);
            return;
        }
        polygons[col-firstCol] = cellPolygons.getTrimmedList();
    }

    public void print() { System.out.println(toString()); }

    public String toString()
    {
        return new String("quadStrip row: " + rowNumber + " firstCol: " + firstCol +
                           " lastCol: " + lastCol + " ID: " + ID);
    }

    public Iterator getPolygonIterator()
    {
        if(polygons == null) polygons = new Object[0][0];
        return new PolygonIterator();

    }

    static public Iterator<QuadStrip> emptyIterator()
    {
        ArrayList<QuadStrip> emptyList  = new ArrayList<QuadStrip>();
        return emptyList.iterator();
    }

    class PolygonIterator implements Iterator
    {
        int i = 0;
        int j = 0;
        int nRows;
        int nCols;
        Object nextPolygon = null;

        PolygonIterator()
        {
            if(polygons == null) polygons = new Object[0][0];
            nRows = polygons.length;
            if(nRows == 0) return;
            nCols = polygons[0].length;
        }

        public boolean hasNext()
        {
            if(j >= nCols)
            {
                i++;
                if(i >= nRows) return false;
                nCols = polygons[i].length;
                j = 0;
                return hasNext();
            }
            nextPolygon = polygons[i][j];
            j++;
            return true;
        }

        public Object next() { return nextPolygon; }

        public void remove() { StsException.systemError(this, "remove", "Cannot remove polygon from quadStrip array."); }
    }
}

