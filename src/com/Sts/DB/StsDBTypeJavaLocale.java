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
public class StsDBTypeJavaLocale extends StsDBTypeClass
{
	public StsDBTypeJavaLocale()
	{
		super();
	}

	public StsDBTypeJavaLocale(int index)
	{
		super(index, java.util.Locale.class);
//        super(index);
    }

	protected Object newInstance()
	{
		return java.util.Locale.getDefault();
	}
}
