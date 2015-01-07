package com.Sts.DB;

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

public class StsDBTypeJavaColor extends StsDBTypeClass
{
	public StsDBTypeJavaColor()
	{
		super();
	}

	public StsDBTypeJavaColor(int index)
	{
		super(index, java.awt.Color.class);
//        super(index);
    }

	protected Object newInstance()
	{
		return new java.awt.Color(0);
	}
}
