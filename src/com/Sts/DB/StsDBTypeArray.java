package com.Sts.DB;

import com.Sts.MVC.*;

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

public class StsDBTypeArray extends StsDBTypeClass implements StsSerializable
{
	private transient int elementTypeIndex;
	private transient byte dimension = 0;
	private transient StsDBType elementType = null;

	private int arrayTypeIndex = -1;

	public StsDBTypeArray()
	{
		super();
	}

	private StsDBTypeArray(int index, Class classType, StsDBObjectTypeList list)
	{
		super(index, classType);
		this.arrayTypeIndex = index;
		initialize(list);
	}

	public StsDBTypeArray(int index)
	{
		super();
		this.arrayTypeIndex = index;
		initializeObjectIO();
	}

	public static StsDBTypeArray getArrayTypeInstance(int index, Class classType, StsDBObjectTypeList list)
	{
		return new StsDBTypeArray(index, classType, list);
	}

	public static StsDBTypeArray getArrayTypeInstance(int index)
	{
		return new StsDBTypeArray(index);
	}

	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		out.writeInt(elementTypeIndex);
		out.writeByte(dimension);
		super.writeObject(out, obj);
	}

	public Object readObject(StsDBInputStream in, Object object) throws IOException
	{
		elementTypeIndex = in.readInt();
		elementType = in.getInputDBType(elementTypeIndex);
		dimension = in.readByte();
		try
		{
			String typeName = Array.newInstance(elementType.getClassType(), new int[dimension]).getClass().getName();
			this.setClassFromTypeName(typeName);
		}
		catch (ClassNotFoundException ex)
		{
		}
		return super.readObject(in, object);
	}

	public Object readObject(StsDBInputStream in) throws IOException
	{
		elementTypeIndex = in.readInt();
		elementType = in.getInputDBType(elementTypeIndex);
		dimension = in.readByte();
		try
		{
			String typeName = Array.newInstance(elementType.getClassType(), new int[dimension]).getClass().getName();
			this.setClassFromTypeName(typeName);
		}
		catch (ClassNotFoundException ex)
		{
		}
		return super.readObject(in);
	}

	protected byte getDimension()
	{
		return dimension;
	}

	public String getTypeName()
	{
		String typeName = super.getTypeName();
		if (typeName == null)
		{
			return "Array";
		}
		return typeName;
	}

	protected StsDBType getElementType()
	{
		return elementType;
	}

    public boolean isArray() { return true; }

    protected boolean initialize(StsDBObjectTypeList list)
	{
		dimension = 0;
		Class c = getClassType();
		while (c.isArray())
		{
			c = c.getComponentType();
			dimension++;
		}
		elementType = list.getDBType(c);
		elementTypeIndex = elementType.getIndex();
		return true;
	}

	public int getIndex()
	{
		return arrayTypeIndex;
	}

	protected StsDBObjectIO getObjectIO()
	{
		return new StsDBObjectIOArray(this);
	}

	/**
	 * copyArrayValues
	 *
	 * @param arrayObject Object
	 * @param size int
	 * @return Object
	 * @throws IllegalAccessException
	 */
	public Object copyArrayValues(Object arrayObject) throws IOException, IllegalAccessException
	{
		return elementType.copyArrayValues(arrayObject);
	}

	/**
	 * copyField
	 *
	 * @param oldObject Object
	 * @param newObject Object
	 * @param field Field
	 * @throws IllegalAccessException
	 */
	public void copyField(Object oldObject, Object newObject, StsDBField field) throws IOException, IllegalAccessException
	{
		Object oldFieldObject = field.getField().get(oldObject);
		Object newFieldObject = copyObject(oldFieldObject);
		field.getField().set(newObject, newFieldObject);
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
		return elementType.readArrayValues(in, arrayObject);
	}

	/**
	 * writeArrayValues
	 *
	 * @param out StsDBOutputStream
	 * @param arrayObject Object
	 * @param size int
	 * @throws IOException
	 */
	public void writeArrayValues(StsDBOutputStream out, Object arrayObject) throws IOException, IllegalAccessException
	{
		elementType.writeArrayValues(out, arrayObject);
	}

	public void addToModel(Object object, StsModel model)
	{
		if (! (elementType instanceof StsDBTypeStsClass)) return;
		addArrayToModel(object, model, dimension);
	}

	private void addArrayToModel(Object arrayObject, StsModel model, int arrayDimension)
	{
		int size = Array.getLength(arrayObject);

		if (arrayDimension == 1)
		{
			for(int n = 0; n < size; n++)
			{
				Object object = Array.get(arrayObject, n);
				elementType.addToModel(object, model);
			}
		}
		else
		{
			for (int n = 0; n < size; n++)
			{
				Object subArray = Array.get(arrayObject, n);
				addArrayToModel(subArray, model, arrayDimension-1);
			}
		}
	}
}
