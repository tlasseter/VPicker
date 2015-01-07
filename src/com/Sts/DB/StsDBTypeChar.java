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

public class StsDBTypeChar extends StsDBType implements StsSerializable
{
	public StsDBTypeChar()
	{
		super(Character.TYPE);
	}

	public StsDBTypeChar(int index)
	{
		super(index, Character.TYPE);
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
		return (char[])arrayObject;
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
		char value = field.getChar(oldObject);
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setChar(newObject, value);
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
		char[] charArray = (char[])arrayObject;
		int size = charArray.length;
		for(int n = 0; n < size; n++)
		{
			charArray[n] = in.readChar();
		}
		return charArray;
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
		char value = in.readChar();
		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.setChar(classObject, value);
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
		char[] charArray = (char[])arrayObject;
		int size = charArray.length;
		for(int n = 0; n < size; n++)
		{
			out.writeChar(charArray[n]);
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
		char value = field.getChar(classObject);
		out.writeChar(value);
	}
}
