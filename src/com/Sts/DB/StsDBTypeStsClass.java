package com.Sts.DB;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

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

public class StsDBTypeStsClass extends StsDBTypeStsSerializable
{
	/** Creates instance of this class with default constructor */
	private transient Constructor constructor = null;

	public StsDBTypeStsClass()
	{
		super();
	}

	public StsDBTypeStsClass(int index, Class classType)
	{
		super(index, classType);
	}

	/** make a copy of this object, or if this is an StsObject and cloneStsObject == true, make a clone (same object, same StsObject.index)  */
	public Object copyObject(Object oldObject)
	{
		try
		{
			StsObject copiedObject = (StsObject) StsObjectCopier.copiedStsObjects.get(oldObject);
			if (copiedObject != null)
			{
				if (debug)
				{
					System.out.println("StsDBTypeStsClass.copyObject() returned existing StsObject " + getTypeName() + "[" + copiedObject.getIndex() + "]");
				}
				return copiedObject;
			}
			StsObject newObject = (StsObject)newInstance();
			if (newObject == null) return null;
			if (debug)
			{
				System.out.println("StsDBTypeStsClass.copyObject() created new StsObject " + getTypeName() + "[" + newObject.getIndex() + "]" +
								   " from " + getTypeName() + "[" + ( (StsObject) oldObject).getIndex() + "]");
			}
			// copyFields will overwrite the newIndex with original index, so restore it
			StsObjectCopier.copiedStsObjects.put(oldObject, newObject);
			boolean skipIndexField = true;
			copyFields(oldObject, newObject, skipIndexField);
			return newObject;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "copyObject(in, obj)", "failed", e);
			return null;
		}
	}

	protected Object newInstance()
	{
		if (constructor == null)
		{
			try
			{
				Class classType = Class.forName(this.getTypeName());
				constructor = classType.getDeclaredConstructor(new Class[] { Boolean.TYPE });
			}
			catch(Exception e)
			{
				StsException.outputWarningException(this, "newInstance", "No default constuctor found for class " + this.getTypeName(), e);
			}
		}
		if (constructor == null)
		{
			return null;
		}
		try
		{
			Object object = constructor.newInstance(new Object[] { Boolean.FALSE });
			if (debug)
			{
				sanityCheck(object);
			}
			return object;
		}
		catch(Exception e)
		{
			StsException.outputException("DBClassType.newInstance() failed for class " + this.getTypeName(), e, StsException.WARNING);
			return null;
		}
	}

	public void addToModel(Object object, StsModel model)
	{
		StsObject stsObject = null;
		try
		{
			// add this object to model
			stsObject = (StsObject)object;
			if(debug) System.out.println("Adding object to model: " + stsObject.getClassname());
			if(stsObject.isPersistent()) return;
			stsObject.addToModel();
			if(debug) System.out.println("    added: " + stsObject.getClassAndNameString());
			stsObject.initialize(model);
			// now add field objects to model
			super.addToModel(object, model);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "addToModel()", "Tried adding object to model: " + stsObject.getClassAndNameString(), e);
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
		Object obj = field.get(classObject);
		if (obj == null)
		{
			out.writeInt(NULL_REFERENCE);
			return;
		}
		out.writeInt(getIndex());

		boolean writeToDisk = field.getDeclaringClass() == StsModel.class; // && obj.getClass() == StsProject.class;
		if (writeToDisk)
		{
			writeObjectFully(out, obj);
		}
		else
		{
			writeObject(out, obj);
		}
	}

	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		int index = ((StsObject)obj).getIndex();
		if (index == -1 ) // && obj.getClass() != StsProject.class)
		{
			System.out.println("Warning, writing StsObject to database with an index of -1. Object = " + obj);
		}
		out.writeInt(index);
	}

	public void writeObjectFully(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		int index = ((StsObject)obj).getIndex();
		if (index == -1) //  && obj.getClass() != StsProject.class)
		{
			System.out.println("Warning, writing StsObject to database with an index of -1. Object = " + obj);
		}
		super.writeObject(out, obj);
	}

	public void exportField(StsDBFileObjectTrader objectTrader, Object classObject, Field field) throws IllegalAccessException
	{
		Object value = field.get(classObject);
		if (value == null)
		{
			return;
		}

		exportObject(objectTrader, value);
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
		// only recurse if object is being added for the first time
		if (objectTrader.add((StsObject)obj))
		{
			if (debug) System.out.println("Adding " + obj + " to list of objects for export");
			super.exportObject(objectTrader, obj);
		}
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
		int classListIndex = in.readInt();
		if (classListIndex == NULL_REFERENCE)
		{
			// A null field indicates that the current definition no longer has that field
			if (field != null)
				field.set(classObject, null);
			return;
		}

		StsDBTypeObject dbClassType =  (StsDBTypeObject)in.getInputDBType(classListIndex);
		if (dbClassType == null)
		{
			throw new NullPointerException("StsDBTypeClass::readField - dbClassType cannot be null");
		}

		Object object = null;
		if (field != null)
		{
			object = field.get(classObject);
		}

		if (field != null && field.getDeclaringClass() == StsModel.class) // && field.getType() == StsProject.class)
		{
			object = dbClassType.readObjectFully(in, object);
		}
		else
		{
			int index = in.readInt();
			Object existingObject = in.getModelObjectOrNull(dbClassType.getClassType(), index);
			if (object == null)
			{
				if (existingObject == null)
				{
					object = in.getModelObject(dbClassType.getClassType(), index);
				}
				else
				{
					object = existingObject;
				}
			}
			else
			{
				if (existingObject == null)
				{
					((StsObject)object).setIndex(index);
					StsObject.getCurrentModel().add((StsObject)object);
				}
				else
				{
					if (existingObject != object)
					{
						if(debug) StsException.systemError(this, "readField", "Existing field object disagrees with object read.");
						object = existingObject;
					}
				}
			}
		}

		// A null field indicates that the current definition no longer has that field
		if (field != null)
			field.set(classObject, object);
	}

	public Object readObjectFully(StsDBInputStream in) throws IOException
	{
		return super.readObject(in);
	}

	public Object readObjectFully(StsDBInputStream in, Object object) throws IOException
	{
		return super.readObject(in, object);
	}

	public Object readObject(StsDBInputStream in) throws IOException
	{
		int index = in.readInt();
		return in.getModelObject(getClassType(), index);
	}

	public Object readObject(StsDBInputStream in, Object object) throws IOException
	{
		int index = in.readInt();
		Object existingObject =  in.getModelObject(getClassType(), index);
		if (existingObject != object)
		{
			if(debug) StsException.systemError(this, "readObject", "Existing object disagrees with object read.");
		}
		return existingObject;
	}
}
