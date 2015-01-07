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

public class StsDBObjectIOByte implements StsDBObjectIO
{
	public StsDBObjectIOByte()
	{
		super();
	}

	/**
	 * copyObject
	 *
	 * @param oldObject Object
	 * @return Object
	 * @throws IllegalAccessException
	 */
	public Object copyObject(Object oldObject) throws IllegalAccessException
	{
		return new Byte(((Byte)oldObject).byteValue());
	}

	/**
	 * readObject
	 *
	 * @param in StsDBInputStream
	 * @return Object
	 * @throws IOException
	 */
	public Object readObject(StsDBInputStream in) throws IOException
	{
		boolean isNull = in.readBoolean();
		if (isNull)
		{
			return null;
		}
		return new Byte(in.readByte());
	}

	/**
	 * readObject
	 *
	 * @param in StsDBInputStream
	 * @param obj Object
	 * @return Object
	 * @throws IOException
	 */
	public Object readObject(StsDBInputStream in, Object obj) throws IOException
	{
		return readObject(in);
	}

	/**
	 * writeObject
	 *
	 * @param out StsDBOutputStream
	 * @param obj Object
	 * @throws IOException
	 * @throws IllegalAccessException
	 */
	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		out.writeBoolean(obj == null);
		if (obj == null)
		{
			return;
		}
		out.writeByte(((Byte)obj).byteValue());
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
	}
}
