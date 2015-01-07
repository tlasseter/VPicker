package com.Sts.Types;

import com.Sts.Utilities.*;

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
public class StsSEGYFormatRec implements Serializable
{
	public String name = null;
	public String userName = null;
	public int loc = 0;
	public byte format = StsSEGYFormat.INT2;
	public String description = null;
	public boolean required = false;
    public boolean userRequired = false;

    public String applyScalar = StsSEGYFormat.scalarAttributes[0];

	public StsSEGYFormatRec()
	{
	}

    public StsSEGYFormatRec(String name, int loc, byte format, String description, boolean required, String applyScalar)
	{
        this(name, loc, format, description, required);
        this.applyScalar = applyScalar;
	}

    public StsSEGYFormatRec(String name, int loc, byte format, String description, boolean required)
    {
        this.name = name;
        this.userName = name;
        this.loc = loc;
        this.format = format;
        this.description = description;
        this.required = required;
	}

	public StsSEGYFormatRec(StsSEGYFormatRec rec)
	{
		this.name = rec.name;
        this.userName = rec.userName;
		this.loc = rec.loc;
		this.format = rec.format;
		this.description = rec.description;
		this.required = rec.required;
        this.applyScalar = rec.applyScalar;
	}

	public void setName(String name)
	{
		this.name = name;
	}

	public String getName()
	{
		return name;
	}

	public boolean setLoc(int loc)
	{
        if(this.loc == loc) return false;
        this.loc = loc;
        return true;
    }

	public int getLoc()
	{
		return loc;
	}

    public void setApplyScalar(String applyScalar)
    {
        this.applyScalar = applyScalar;
    }

    public boolean setApplyScalarSelectedIndex(int selectedIndex)
    {
        if(this.applyScalar == StsSEGYFormat.scalarAttributes[selectedIndex]) return false;
        this.applyScalar = StsSEGYFormat.scalarAttributes[selectedIndex];
        return true;
    }

    public String getApplyScalar()
    {
        return applyScalar;
	}
    public boolean getRequired()
    {
        return required;
    }

    public void setRequired(boolean required)
    {
        this.required = required;
    }

    public void setUserName(String name)
    {
        this.userName = name;
    }

    public void setDescription(String desc)
    {
        this.description = desc;
    }

    public String getUserName()
    {
        return userName;
    }

    public String getDescription()
    {
        return description;
    }

	public Object clone() throws CloneNotSupportedException
	{
		StsSEGYFormatRec tmp = (StsSEGYFormatRec) super.clone();
		return tmp;
	}
/*
	public void setHdrValue(byte[] hdr, double val)
	{
		formatChanged = true;
		setHdrValue(hdr, val, littleEndian);
		return;
	}
*/
/*
	public void setHdrValue(byte[] hdr, double val, boolean littleEndian)
	{
		int tInt = 0;
		short sInt = 0;

		switch(format)
		{
		   case StsSEGYFormat.INT2 : // 2 byte integer
			   tInt = (short) val;
//                   hdr[loc] =
			   break;
		   case StsSEGYFormat.INT4 : // 4 byte integer
			   tInt = (int) val;
//                   hdr[loc] =
			   break;
			default:
				;
		}
	}
*/
	// return a double gives us accuracy to cover complete range: byte to float
	public double getHdrValue(StsSegyVolume segyVolume)
	{
		byte[] hdr = segyVolume.getBinaryHeader();
		boolean isLittleEndian = segyVolume.getIsLittleEndian();
		return getHdrValue(hdr, isLittleEndian);
	}

	// return a double gives us accuracy to cover complete range: byte to float
	public double getHdrValue(byte[] hdr, boolean littleEndian)
	{
		switch(format)
		{
		   case StsSEGYFormat.BYTE : // 8 bit integer
				return (double)hdr[loc];
		   case StsSEGYFormat.INT2 : // 2 byte integer
				return (double) StsMath.convertBytesToShort(hdr, loc, littleEndian);
		   case StsSEGYFormat.INT4 : // 4 byte integer
				return (double)StsMath.convertIntBytes(hdr, loc, littleEndian);
		   case StsSEGYFormat.IBMFLT : // 4 byte IBM floating point
				return (double)StsMath.convertIBMFloatBytes(hdr, loc, littleEndian);
		   case StsSEGYFormat.IEEEFLT : // 4 byte IEEE Float
				return (double)Float.intBitsToFloat(StsMath.convertIntBytes(hdr, loc, littleEndian));
			default:
				return -1;
		}
	}

	public double getHdrValue(byte[] hdr, int offset, boolean littleEndian)
	{
		switch(format)
		{
		   case StsSEGYFormat.BYTE : // 8 bit integer
				return (double)hdr[loc + offset];
		   case StsSEGYFormat.INT2 : // 2 byte integer
				return (double) StsMath.convertBytesToShort(hdr, loc + offset, littleEndian);
		   case StsSEGYFormat.INT4 : // 4 byte integer
				return (double)StsMath.convertIntBytes(hdr, loc + offset, littleEndian);
		   case StsSEGYFormat.IBMFLT : // 4 byte IBM floating point
				return (double)StsMath.convertIBMFloatBytes(hdr, loc + offset, littleEndian);
		   case StsSEGYFormat.IEEEFLT : // 4 byte IEEE Float
				return (double)Float.intBitsToFloat(StsMath.convertIntBytes(hdr, loc + offset, littleEndian));
			default:
				return -1;
		}
	}

    public double getHdrValue(byte[] hdr, int offset, boolean littleEndian, int xyScalar, int edScale)
    {
        double scale = 1.0f;
        if(applyScalar.equals("CO-SCAL"))
        {
            if (xyScalar < 0)
                scale = -(double)1/(double)xyScalar;
            else if (xyScalar > 0)
                scale = (double)xyScalar;
        }
        else if(applyScalar.equals("ED-SCAL"))
        {
            if (edScale < 0)
               scale = -(double)1/(double)edScale;
           else if (edScale > 0)
               scale = (double)edScale;
        }
        switch(format)
        {
           case StsSEGYFormat.BYTE : // 8 bit integer
                return (double)hdr[loc + offset] * scale;
           case StsSEGYFormat.INT2 : // 2 byte integer
                return (double) StsMath.convertBytesToShort(hdr, loc + offset, littleEndian) * scale;
           case StsSEGYFormat.INT4 : // 4 byte integer
                return (double)StsMath.convertIntBytes(hdr, loc + offset, littleEndian) * scale;
           case StsSEGYFormat.IBMFLT : // 4 byte IBM floating point
                return (double)StsMath.convertIBMFloatBytes(hdr, loc + offset, littleEndian) * scale;
           case StsSEGYFormat.IEEEFLT : // 4 byte IEEE Float
                return (double)Float.intBitsToFloat(StsMath.convertIntBytes(hdr, loc + offset, littleEndian)) * scale;
            default:
                return -1;
        }
	}

	public boolean setFormat(int format)
	{
		if (this.format == (byte) format)return false;
		this.format = (byte) format;
		return true;
	}

	public int getFormat()
	{
		return format;
	}

    public boolean isRequired()
    {
        return required || userRequired;
    }
}
