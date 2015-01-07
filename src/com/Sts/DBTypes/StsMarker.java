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

public class StsMarker extends StsMainObject implements StsSelectable
{
	// constants
	static public final byte GENERAL = 0;
	static public final byte SURFACE = 1;
	static public final byte FAULT = 2;
	static public final byte UNCONFORMITY = 3;
	static public final byte EQUIPMENT = 4;
	static public final byte PERFORATION = 5;
	static public final byte FMI = 6;
    static public final byte OFFSET_SURFACE = 7;

    static public final String[] markerTypes = new String[] { "General", "Surface", "Fault", "Unconformity", "Equipment", "Perforation", "FMI", "OffsetSurface"};
	
    // instance fields
	protected StsColor stsColor = null;
	protected StsSurface surface = null;
    protected StsModelSurface modelSurface = null;
	protected StsObjectRefList wellMarkers = null;
	protected StsWellZoneSet wellZoneSetAbove = null;
	protected StsWellZoneSet wellZoneSetBelow = null;

    static int colorIndex = 0;

	/** constructor */
	public StsMarker()
	{
	}

	public StsMarker(boolean persistent)
	{
		super(persistent);
	}

	public StsMarker(String name, byte type)
	{
        super(false);
        this.type = type;
        setName(name);
		initializeColor();
        addToModel();
	}

    public StsMarker(String name, byte type, StsSurface surface)
    {
        super(false);
        this.type = type;
        setName(name);
        if(surface != null)
			setSurface(surface);
		else
		    initializeColor();
        addToModel();
    }

    public StsMarker(String name, byte type, StsColor color)
    {
        super(false);
        this.type = type;
        setName(name);
        setStsColor(color);
        addToModel();
    }
    
	private void initializeColor()
	{
		StsSpectrum spectrum = currentModel.getSpectrum("Basic");
		StsColor color = spectrum.getColor(colorIndex++);
        setStsColor(color);
	}

    public boolean initialize(StsModel model)
    {
        return true;
    }

	/** accessors */

	public String getName()
	{return super.getName();
	}

	public void setStsColor(StsColor color)
	{
		if(stsColor != null && stsColor.equals(color)) return;
		stsColor = new StsColor(color);
		if(isPersistent()) dbFieldChanged("stsColor", color);
	}

	public StsColor getStsColor()
	{
		if(surface != null) return surface.getStsColor();
		return stsColor;
	}

    public void setSurface(StsSurface surface)
	{
		this.surface = surface;
		if (surface != null) setStsColor(surface.getStsColor());
        if(isPersistent()) fieldChanged("surface", surface);
	}

	public StsSurface getSurface()
	{
        return surface;
	}

    public void setModelSurface(StsModelSurface surface)
    {
        this.modelSurface = surface;
        this.type = SURFACE;
        dbFieldChanged("type", SURFACE);
        dbFieldChanged("modelSurface", modelSurface);
        if (surface != null) setStsColor(surface.getStsColor());
    }

    public StsModelSurface getModelSurface()
    {
        return modelSurface;
    }

	public void setWellZoneSetAbove(StsWellZoneSet wellZoneSet)
	{
		wellZoneSetAbove = wellZoneSet;
	}

	public StsWellZoneSet getWellZoneSetAbove()
	{return wellZoneSetAbove;
	}

	public void setWellZoneSetBelow(StsWellZoneSet wellZoneSet)
	{
		wellZoneSetBelow = wellZoneSet;
	}

	public StsWellZoneSet getWellZoneSetBelow()
	{return wellZoneSetBelow;
	}

	/** check type for valid value */
	static public boolean isTypeValid(byte type)
	{
        return type >= 0 && type < 8;
	}

	/** get a string equivalent for a type */
	static public String typeToString(byte type)
	{
        if(type < 0 || type >= 8) return null;
        return markerTypes[type];
	}

    /** get a string equivalent for a type */
    static public byte stringToType(String typeString)
    {
        for(int n = 0; n < 8; n++)
            if(typeString.equals(markerTypes[n])) return (byte)n;
        return (byte)-1;
    }

	/** add a well marker to the reflist */
	public void addWellMarker(StsWellMarker wellMarker)
	{
		if (wellMarkers == null)
		{
			wellMarkers = StsObjectRefList.constructor(1, 1, "wellMarkers", this);
		}
		wellMarkers.add(wellMarker);
	}
	
	/** add a well marker to the reflist */
/*	public void addPerforationMarker(StsPerforationMarker perfMarker)
	{
		if (wellMarkers == null)
		{
			wellMarkers = StsObjectRefList.constructor(1, 1, "wellMarkers", this);
		}
		wellMarkers.add(perfMarker);
	}
*/	
	public StsObjectRefList getWellMarkers()
	{
		return wellMarkers;
	}
	
	public void removeWellMarkers(StsWellMarker wellMarker) throws StsException
	{
		if (wellMarker == null)
		{
			return;
		}
		wellMarkers.delete(wellMarker);
	}

	/** get average marker z value to use in ordering markers */
	public float getOrderingValue()
	{
		try
		{
			int nWellMarkers = wellMarkers.getSize();
			if (nWellMarkers < 1)
			{
				return StsParameters.nullValue;
			}
			int nValues = 0;
			float totalValues = 0.0f;
			for (int i = 0; i < nWellMarkers; i++)
			{
				StsWellMarker wellMarker = (StsWellMarker) wellMarkers.getElement(i);
				float value = wellMarker.getZ();
				if (value == StsParameters.nullValue)
				{
					continue;
				}
				nValues++;
				totalValues += value;
			}
			if (nValues == 0)
			{
				return StsParameters.nullValue;
			}
			return totalValues / (float) nValues;
		}
		catch (Exception e)
		{return StsParameters.nullValue;
		}
	}

	static public void tryToLinkAdjacentIntervals(StsMarker marker, StsModelSurface surface)
	{
		tryToLinkIntervalAbove(marker, surface);
		tryToLinkIntervalBelow(marker, surface);
	}

	static public void tryToLinkIntervalAbove(StsMarker marker, StsModelSurface surface)
	{
		if (!areLinked(marker, surface))
		{
			return;
		}

		StsWellZoneSet wellZoneAbove = marker.getWellZoneSetAbove();
		StsZone zoneAbove = surface.getZoneAbove();
		if (!StsWellZoneSet.areLinked(wellZoneAbove, zoneAbove))
		{
			StsMarker topMarker = (wellZoneAbove == null) ? null : wellZoneAbove.getTopMarker();
			StsModelSurface topSurface = (zoneAbove == null) ? null : zoneAbove.getTopModelSurface();
			if (areLinked(topMarker, topSurface))
			{
				zoneAbove.setWellZoneSet(wellZoneAbove);
				wellZoneAbove.setParentZone(zoneAbove);
			}
		}
	}

	static public void tryToLinkIntervalBelow(StsMarker marker, StsModelSurface surface)
	{
		if (!areLinked(marker, surface))
		{
			return;
		}

		StsWellZoneSet wellZoneBelow = marker.getWellZoneSetBelow();
		StsZone zoneBelow = surface.getZoneBelow();
		if (!StsWellZoneSet.areLinked(wellZoneBelow, zoneBelow))
		{
			StsMarker baseMarker = (wellZoneBelow == null) ? null : wellZoneBelow.getBaseMarker();
			StsModelSurface botSurface = (zoneBelow == null) ? null : zoneBelow.getBaseModelSurface();
			if (areLinked(baseMarker, botSurface))
			{
				zoneBelow.setWellZoneSet(wellZoneBelow);
				wellZoneBelow.setParentZone(zoneBelow);
			}
		}
	}

	static public boolean areLinked(StsMarker m, StsModelSurface s)
	{
		if (m == null || s == null)
		{
			return false;
		}
		if (s.getMarker() != m)
		{
			return false;
		}
		return true;
	}

	static public void tryToBuildNewIntervals(StsModel model, StsMarker m, StsModelSurface s)
	{
		// see if we have a null horizon
		if (m == null || s == null)
		{
			return;
		}

		// be sure the marker & surface are correlated
		if (m.getSurface() != s || s.getMarker() != m)
		{
			return;
		}

		StsWellZoneSet wzAbove = m.getWellZoneSetAbove();
		StsWellZoneSet wzBelow = m.getWellZoneSetBelow();
		StsZone zAbove = s.getZoneAbove();
		StsZone zBelow = s.getZoneBelow();
		try
		{
			if (wzAbove == null && zAbove != null)
			{
				StsMarker mAbove = zAbove.getTopModelSurface().getMarker();
				if (mAbove != null)
				{
					wzAbove = new StsWellZoneSet(StsWellZoneSet.STRAT, mAbove, m);
					wzAbove.setParentZone(zAbove);
					zAbove.setWellZoneSet(wzAbove);
					StsWellZoneSet.buildWellZones(model, mAbove, m);
				}
			}
			else if (wzAbove != null && zAbove == null)
			{
				StsModelSurface sAbove = wzAbove.getTopMarker().getModelSurface();
				if (sAbove != null)
				{
					zAbove = new StsZone(wzAbove);
					zAbove.setTopModelSurface(sAbove);
					zAbove.setBaseModelSurface(m.getModelSurface());
				}
			}

			if (wzBelow == null && zBelow != null)
			{
				StsMarker mBelow = zBelow.getBaseModelSurface().getMarker();
				if (mBelow != null)
				{
					wzBelow = new StsWellZoneSet(StsWellZoneSet.STRAT, m, mBelow);
					wzBelow.setParentZone(zBelow);
					zBelow.setWellZoneSet(wzBelow);
					StsWellZoneSet.buildWellZones(model, m, mBelow);
				}
			}
			else if (wzBelow != null && zBelow == null)
			{
				StsModelSurface sBelow = wzBelow.getBaseMarker().getModelSurface();
				if (sBelow != null)
				{
					zBelow = new StsZone(wzBelow);
					zBelow.setTopModelSurface(m.getModelSurface());
					zBelow.setBaseModelSurface(sBelow);
				}
			}
		}
		catch (Exception e)
		{}
	}
}