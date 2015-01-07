package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Types.*;
/**
 * Stores an OpenSpirit-derived formation marker
 * @author lkw
 *
 */
public class StsOSWellMarker
{
    // constants
    static public final byte UNSPECIFIED_TYPE = 0;
    static public final byte TOP = 1;
    static public final byte BASE = 2;
    static public final byte FAULT = 3;
    static public final byte UNCONFORMITY = 4;
    
    // attributes
    private String name = null;
    private byte type = UNSPECIFIED_TYPE;
    private float mdepth = 0.f;
    private float tvdepth = 0.f;
    private float twt = 0.f;
    private StsPoint location = null;
    
    public StsOSWellMarker(String name, byte type, float mdepth, float tvdepth, 
    		float twt, double x, double y)
    {
    	this.name = name;
    	this.type = type;
    	this.mdepth = mdepth;
    	this.tvdepth = tvdepth;
    	this.twt = twt;
    	location = new StsPoint(x, y, tvdepth);
    }

	public StsPoint getLocation()
	{
		return location;
	}

	public void setLocation(StsPoint location)
	{
		this.location = location;
	}

	public float getMdepth()
	{
		return mdepth;
	}

	public void setMdepth(float mdepth)
	{
		this.mdepth = mdepth;
	}

	public String getName()
	{
		return name;
	}

	public void setName(String name)
	{
		this.name = name;
	}

	public float getTvdepth()
	{
		return tvdepth;
	}

	public void setTvdepth(float tvdepth)
	{
		this.tvdepth = tvdepth;
	}

	public float getTwt()
	{
		return twt;
	}

	public void setTwt(float twt)
	{
		this.twt = twt;
	}

	public byte getType()
	{
		return type;
	}

	public void setType(byte type)
	{
		this.type = type;
	}
}
