package com.Sts.DB;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsDBTypeJavaColorSpace extends StsDBTypeClass
{
	public StsDBTypeJavaColorSpace()
	{
		super();
	}

	public StsDBTypeJavaColorSpace(int index)
	{
//		super(index);
        super(index, java.awt.color.ColorSpace.class);
	}

	protected Object newInstance()
	{
		return java.awt.color.ColorSpace.getInstance(java.awt.color.ColorSpace.CS_LINEAR_RGB);
	}

    public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		super.writeObject(out, obj);
	}

    public Object readObject(StsDBInputStream in, Object object) throws IOException
	{
        return super.readObject(in, object);
    }
}
