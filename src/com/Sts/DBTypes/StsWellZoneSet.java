
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.text.*;

public class StsWellZoneSet extends StsMainObject implements StsSelectable
{
    // constants
    static public final int STRAT = 0;
    static public final int LITH = 1;

    // instance fields
    protected int zoneType = STRAT;
    protected int nSubZones = StsZone.DEFAULT_N_SUBZONES;
    protected int subZoneType = StsZone.SUBZONE_UNIFORM;
    protected StsColor stsColor = null;
    protected StsObjectRefList wellZones;
    protected boolean displayOnWellLines = true;
    protected StsMarker topMarker = null;
    protected StsMarker baseMarker = null;
    protected StsZone parentZone;

    /** constructor */
    public StsWellZoneSet(int type, String name)
    {
        if (name==null)
        {
            StsException.systemError("StsWellZoneSet.StsWellZoneSet:"
                      + " Cannot create a well zone set with a null name.");
        }
        setName(name);
        setZoneType(type);
    }

    /** constructor */
    public StsWellZoneSet(int type, StsMarker topMarker, StsMarker baseMarker)
            throws StsException
    {
        if (topMarker==null || baseMarker==null)
        {
            StsException.systemError("StsWellZone.StsWellZone:"
                      + " Cannot create a well zone with a null top or base marker.");
        }
        setZoneType(type);
        setName(topMarker.getName());
        stsColor = topMarker.getStsColor();
        this.topMarker = topMarker;
        topMarker.setWellZoneSetBelow(this);
        this.baseMarker = baseMarker;
        baseMarker.setWellZoneSetAbove(this);
    }

    /** DB constructor */
	public StsWellZoneSet()
	{
	}

    public boolean initialize(StsModel model) { return true; }

    // Accessors
    public boolean setZoneType(int type)
    {
        switch(type)
        {
            case STRAT:
            case LITH:
                this.zoneType = type;
                return true;
        }
        return false;
    }
    public int getZoneType(){ return zoneType; }
    public void setStsColor(StsColor color) { setStsColor(color, true); }
    public void setStsColor(StsColor color, boolean setParent)
    {
        if (stsColor.equals(color)) return;  // no change
            stsColor = color;
        if (parentZone!=null && setParent)
            parentZone.setStsColor(color);

        // need to rebuild all well lines with new color
        setWellLinesNeedRebuild();
    }
    public StsColor getStsColor() { return stsColor; }

    public void setTopMarker(StsMarker marker) { topMarker = marker; }
    public StsMarker getTopMarker() { return topMarker; }
    public void setBaseMarker(StsMarker marker) { baseMarker = marker; }
    public StsMarker getBaseMarker() { return baseMarker; }

    public void setParentZone(StsZone parentZone)
    {
        this.parentZone = parentZone;
        if (parentZone!=null)
        {
            // set data from parent zone
            setStsColor(parentZone.getStsColor(), false);
            setNSubZones(parentZone.getNSubZones(), false);
            setSubZoneType(parentZone.getSubZoneType(), false);
        }
    }
    public StsZone getParentZone() { return parentZone; }

    /** sub zone methods */
    public boolean setNSubZones(int nSubZones, boolean setParent)
    {
        if (nSubZones <= 0) return false;
        this.nSubZones = nSubZones;

        // set zone if requested
        if (parentZone!=null && setParent)
        {
//           if (!parentZone.setNSubZones(nSubZones)) return false;
        }
        return true;
    }
    public int getNSubZones() { return nSubZones; }

    public boolean setSubZoneType(int type, boolean setParent)
    {
        switch(type)
        {
            case StsZone.SUBZONE_UNIFORM:
            case StsZone.SUBZONE_ONLAP:
            case StsZone.SUBZONE_OFFLAP:
            case StsZone.SUBZONE_TRUNCATED:
                subZoneType = type;
                // set zone if requested
                if (parentZone == null || !setParent) return true;
                return parentZone.setSubZoneType(type);
            default:
                return false;
        }
    }

    public int getSubZoneType() { return subZoneType; }

    /** add a well zone to the reflist */
    public void addWellZone(StsWellZone wellZone)
    {
        if (wellZone==null) return;
        if(wellZones == null)
        {
            wellZones = StsObjectRefList.constructor(1, 1, "wellZones", this);
        }
        wellZones.add(wellZone);
    }
    public StsObjectRefList getWellZones() { return wellZones; }
    public void removeWellZone(StsWellZone wellZone) throws StsException
    {
        if (wellZone==null) return;
        wellZones.delete(wellZone);
    }

    /** display subzone properties for all the well zones & a log name */
    public void displayWellSubZoneProperties(String logName)
    {
        StringBuffer buffer = new StringBuffer("Well Zone Set: " + getName() + "\n \n");
        buffer.append("Well sub-zone properties (bottom to top)\n");
        buffer.append("Log: " + logName + " \n \n");

        if (logName==null)
        {
            buffer.append("Null log curve name.\n");
        }
        else if (wellZones==null)
        {
            buffer.append("No well zones found!\n");
        }
        else
        {
            NumberFormat fmt1 = NumberFormat.getNumberInstance();
            fmt1.setMaximumFractionDigits(2);
            fmt1.setMinimumFractionDigits(2);
            fmt1.setGroupingUsed(false);
            NumberFormat fmt2 = NumberFormat.getNumberInstance();
            fmt2.setMaximumFractionDigits(5);
            fmt2.setMinimumFractionDigits(5);
            fmt2.setGroupingUsed(false);
            int nWellZones = wellZones.getSize();
            boolean averagesFound = false;
            for (int i=0; i<nWellZones; i++)
            {
                int index = 0;
                StsWellZone wellZone = (StsWellZone)wellZones.getElement(i);
                float[] subZoneAverages = wellZone.getSubZoneAverages(logName);
                if (subZoneAverages==null) continue;
                averagesFound = true;
                buffer.append("Well: " + wellZone.getWell().getName() + " \n");
                float[] subZoneCenterZs = wellZone.getLayerCenterZs();
                float[] subZoneThicknesses = wellZone.getLayerThicknesses();
                // go from base to top
                buffer.append("\tindex\tdepth\t\tthick\taverage\n");
                for (int j=subZoneCenterZs.length-1; j>=0; j--)
                {
                    StringBuffer centerZ = new StringBuffer(fmt1.format(subZoneCenterZs[j]));
                    if (centerZ.length()==7) centerZ.insert(0, ' ');
                    StringBuffer thick = new StringBuffer(fmt1.format(subZoneThicknesses[j]));
                    if (thick.length()==4) thick.insert(0, ' ');
                    buffer.append("\t" + index + "\t" + centerZ.toString()
                            + "\t" + thick.toString()
                            + "\t" + fmt2.format(subZoneAverages[j]) + "\n");
                    index++;
                }
                buffer.append(" \n");
            }
            if (!averagesFound) buffer.append("No subzone averages found!\n \n");
        }

        // print out buffer
        PrintWriter out = new PrintWriter(System.out, true); // needed for correct formatting
        out.println(buffer.toString());

        // display dialog box
        StsTextAreaDialog dialog = new StsTextAreaDialog(null, "Zone properties Listing for "
                + getName(), buffer.toString(), 30, 50, false);
        dialog.setVisible(true);
    }

    /** set/get flag for display of associated well zones on well lines */
    public void setDisplayOnWellLines(boolean state) { displayOnWellLines = state; }
    public boolean getDisplayOnWellLines() { return displayOnWellLines; }
    public void setWellLinesNeedRebuild()
    {
        if (wellZones==null) return;
        int nWellZones = wellZones.getSize();
        for (int i=0; i<nWellZones; i++)
        {
            StsWellZone wellZone = (StsWellZone)wellZones.getElement(i);
            if (wellZone==null) continue;
            StsWell well = wellZone.getWell();
//            if (well!=null) well.setWellLineNeedsRebuild();
        }
    }

    public float getOrderingValue() { return (float)wellZones.getSize(); }

    static public boolean areLinked(StsWellZoneSet wzs, StsZone z)
    {
        if (wzs == null || z == null) return false;
        if (z.getWellZoneSet() != wzs) return false;
        return true;
    }

    static public boolean isNewWellZoneSet(StsModel model, StsMarker top, StsMarker base)
    {
        try
        {
            StsClass wellZoneSets = model.getCreateStsClass(StsWellZoneSet.class);
            int nWellZoneSets = (wellZoneSets==null) ? 0 : wellZoneSets.getSize();
            for (int i=0; i<nWellZoneSets; i++)
            {
                StsWellZoneSet wellZoneSet = (StsWellZoneSet)wellZoneSets.getElement(i);
                if (wellZoneSet==null) continue;  // shouldn't happen
                if (top==wellZoneSet.getTopMarker()) return false;    // already used as a top of wellZone
                if (base==wellZoneSet.getBaseMarker()) return false;  // already used as a base of wellZone
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return false;
        }
    }

    public boolean buildZone(StsModel model)
    {
        if (topMarker == null || baseMarker == null) return false;

        StsModelSurface topSurface = topMarker.getModelSurface();
        StsModelSurface baseSurface = baseMarker.getModelSurface();
        if (topSurface == null || baseSurface == null) return false;
        StsZoneClass zoneClass = (StsZoneClass)model.getCreateStsClass(StsZone.class);
        if(zoneClass.getExistingZone(topSurface, baseSurface) == null)
        {
            StsZone zone = null;
            try { zone = new StsZone(topSurface.getName(), topSurface, baseSurface); }
            catch (StsException e) { return false; }
            setParentZone(zone);
            zone.setWellZoneSet(this);
            return true;
        }

        return false;
    }

    public static boolean buildWellZones(StsModel model, StsMarker top,
            StsMarker base)
    {
        try
        {
            StsObjectRefList topWellMarkers = top.getWellMarkers();
            StsObjectRefList baseWellMarkers = base.getWellMarkers();
            if (topWellMarkers==null || baseWellMarkers==null) return false;
            int nTopWellMarkers = topWellMarkers.getSize();
            int nBaseWellMarkers = baseWellMarkers.getSize();
            if (nTopWellMarkers==0 || nBaseWellMarkers==0) return false;
            for (int i=0; i<nTopWellMarkers; i++)
            {
                StsWellMarker topWellMarker = (StsWellMarker)topWellMarkers.getElement(i);
                if (topWellMarker==null) return false;  // shouldn't happen
                StsWell well = topWellMarker.getWell();
                for (int j=0; j<nBaseWellMarkers; j++)
                {
                    StsWellMarker baseWellMarker = (StsWellMarker)baseWellMarkers.getElement(j);
                    if (well==baseWellMarker.getWell())  // found top & base for same well
                    {
                        try
                        {
                            if (topWellMarker.getZ() >= baseWellMarker.getZ()) break;
                            new StsWellZone(well, StsWellZoneSet.STRAT,
                                    topWellMarker.getName(), topWellMarker,
                                    baseWellMarker, null);
                        }
                        catch (Exception e) { return false; }
                        break;
                    }
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return false;
        }
    }
}













