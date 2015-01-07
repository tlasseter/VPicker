package com.Sts.DB;

import java.io.*;
import java.lang.reflect.*;

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

public class StsDBTypeByte extends StsDBType implements StsSerializable
{
	public StsDBTypeByte()
	{
		super(Byte.TYPE);
	}

	public StsDBTypeByte(int index)
	{
		super(index, Byte.TYPE);
	}

	/**
	 * copyArrayValues
	 *
	 * @param arrayObject Object
	 * @param size int
	 * @return Object
	 * @throws IllegalAccessException
	 */
	public Object copyArrayValues(Object arrayObject) throws IllegalAccessException
	{
		return (byte[])arrayObject;
	}

	/**
	 * copyField
	 *
	 * @param oldObject Object
	 * @param newObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 */
	public void copyField(Object oldObject, Object newObject, Field field) throws IllegalAccessException
	{
		byte value = field.getByte(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setByte(newObject, value);
	}

	/**
	 * readArrayValues
	 *
	 * @param in StsDBInputStream
	 * @param size int
	 * @return Object
	 * @throws IOException
	 */
	public Object readArrayValues(StsDBInputStream in, Object arrayObject) throws IOException
	{
		byte[] byteArray = (byte[])arrayObject;
		int size = byteArray.length;
		for(int n = 0; n < size; n++)
		{
			byteArray[n] = in.readByte();
		}
		return byteArray;
	}

	/**
	 * readField
	 *
	 * @param in StsDBInputStream
	 * @param classObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 * @throws IOException
	 */
	public void readField(StsDBInputStream in, Object classObject, Field field) throws IllegalAccessException, IOException
	{
		byte value = in.readByte();
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setByte(classObject, value);
	}

	/**
	 * writeArrayValues
	 *
	 * @param out StsDBOutputStream
	 * @param arrayObject Object
	 * @param size int
	 * @throws IOException
	 */
	public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException
	{
		byte[] byteArray = (byte[])arrayObject;
		int size = byteArray.length;
		for(int n = 0; n < size; n++)
		{
			out.writeByte(byteArray[n]);
		}
	}

	/**
	 * writeField
	 *
	 * @param out StsDBOutputStream
	 * @param classObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 * @throws IOException
	 */
	public void writeField(StsDBOutputStream out, Object classObject, Field field) throws IllegalAccessException, IOException
	{
		byte value = field.getByte(classObject);
		out.writeByte(value);
	}
}
