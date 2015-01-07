package com.Sts.DB;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.util.*;

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

public class StsObjectCopier
{
	StsModel model;
	StsObject object;
	StsDBTypeStsClass classType;

	/** DBTypes list (primitives, arrays, and classes) used in copy operations */
	static public StsDBObjectTypeList typeList = new StsDBObjectTypeList("objectCopier");
	/** Contains oldObject (key) and newObject (value) for StsObject copy operations.
	 *  Any subsequent copy of oldObject will use this newObject instead.
	 */
	static public HashMap copiedStsObjects = new HashMap();

	static final boolean debug = false;

	public StsObjectCopier(StsModel model, StsObject object)
	{
		this.model = model;
		this.object = object;
		classType = (StsDBTypeStsClass)typeList.getDBType(object);
//		DBClassType dbType = DBClassType.constructor(object.getClass(), typeList);
	}

	static public StsDBObjectTypeList getDBTypeList() { return typeList; }

	public boolean setFieldCloneStsObject(String fieldName)
	{
		return setFieldCloneStsObjectFlag(fieldName, true);
	}

	public boolean setFieldCopyStsObject(String fieldName)
	{
		return setFieldCloneStsObjectFlag(fieldName, false);
	}

	private boolean setFieldCloneStsObjectFlag(String fieldName, boolean clone)
	{

		StsDBField fieldType = classType.getDBFieldType(fieldName);
		if(fieldType == null) return false;
		StsDBType fieldClassType = fieldType.getDBType();
		if(!(fieldClassType instanceof StsDBTypeObject))
		{
			StsException.systemError("StsObjectCopier.setFieldCloneStsObject() failed. " + fieldType.getFieldName() + " is not a class.");
			return false;
		}
		else if(fieldClassType instanceof StsDBTypeArray)
		{
			StsDBType elementType = ((StsDBTypeArray)fieldClassType).getElementType();
			if(!(elementType instanceof StsDBTypeStsClass))
				StsException.systemError("StsObjectCopier.setFieldCloneStsObject() failed. " + fieldType.getFieldName() + " is not an array of StsObjects.");
			return false;
		}
		else if(!(classType instanceof StsDBTypeStsClass))
		{
			StsException.systemError("StsObjectCopier.setFieldCloneStsObject() failed. " + fieldType.getFieldName() + " is not an StsObject class.");
			return false;
		}
		fieldType.setCloneFlag(clone);
		return true;
	}

	public boolean cloneAllFields()
	{
		return classType.cloneAllFields();
	}

	public Object copy()
	{
		StsObject newObject = (StsObject)classType.copyObject(object);
		classType.addToModel(newObject, model);
		return newObject;
	}

	static public void main(String[] args)
	{
		testObjectCopy();
//		testDBWriteRead("c:\\stsdev\\c75n-tom\\dbWriteReadTest");
	}

	static private void testObjectCopy()
	{
		// setup the model for testing
		StsModel model = StsModel.constructor();
		// create an object to copy
		StsObjectTest object = new StsObjectTest();
		object.initialize();
		object.setName("original before copy");
		// output original before copy to compare with after copy
		object.print();
		// copy the object
		StsObjectCopier copier = new StsObjectCopier(model, object);
//		copier.setFieldCloneStsObject("refList");
		StsObjectTest objectCopy = (StsObjectTest)copier.copy();
		// futz with some of the new copy fields
		objectCopy.setName("copy");
		StsObjectFieldTest element = (StsObjectFieldTest)objectCopy.refList.getElement(0);
		if(element != null) element.i = 20;
		// print out both to see how they compare
		object.setName("original after copy");
		object.print();
		objectCopy.print();
		String pathname = "c:\\stsdev\\c80_54\\dbtest";
		StsDBFileObject.writeObjectFile(pathname, object, null);
		StsModel copyModel = StsModel.constructor();
        object = new StsObjectTest();
        StsDBFileObject.readDatabaseObjectFile(pathname, object, null);
        StsObject[] stsObjects = copyModel.getObjectList(StsObjectTest.class);
        for(StsObject stsObject : stsObjects)
            ((StsObjectTest)stsObject).print();    
		System.exit(0);
	}
/*
	static private void testDBWriteRead(String pathname)
	{
		StsModel model = new StsModel();
		model.project = new StsProject();
		model.initializeModel();
		model.setDatabase(new StsDBFile());
		StsObjectTest test = new StsObjectTest();
		test.initializeSuperGather();
		test.setName("dbTest");
		test.print();
		StsDBFileNew.writeObjectFile(pathname, test);
		test = (StsObjectTest)StsDBFileNew.readDatabaseObjectFile(pathname);
		test.print();
	}
*/
}

class StsObjectTest extends StsObject
{
	String name = null;
	StsObjectRefList refList = null;
// 	StsObjectRefList colorscales = null;

	StsObjectTest()
	{
	}

	StsObjectTest(boolean persistent)
	{
		super(persistent);
	}

	public boolean initialize(StsModel model)
	{
		return true;
	}

	void initialize()
	{
		StsObjectRefList refList = StsObjectRefList.constructor(1, 1, "refList", this);
		StsObjectFieldTest objectFieldTest = new StsObjectFieldTest();
		objectFieldTest.setParent(this);
		refList.add(objectFieldTest);
        /*
		StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);
		StsColorscale colorscale = new StsColorscale("FRED", spectrum, 0.0f, 255f);
		if(colorscales == null) colorscales = StsObjectRefList.constructor(2, 2, "colorscales", this);
		colorscales.add(colorscale);
		*/
	}

	void setName(String name)
	{
		this.name = name;
	}

    void print()
	{
		System.out.println("StsObjectTest " + name + " index: " + getIndex());
		if(refList == null)
			System.out.println("    refList is null");
		else
		{
			int nElements = refList.getSize();
			System.out.println("    refList " + refList.getName() + " nElements " + nElements);
			for (int n = 0; n < nElements; n++)
			{
				StsObjectFieldTest element = (StsObjectFieldTest) refList.getElement(n);
                if (element != null)
				{
					System.out.println("        refList element index " + n + " i " + element.i);
					System.out.println("        refList element parent index " + element.parent.getIndex() + " should be " + getIndex());
				}
				else
					System.out.println("        refList element index " + n + " i " + element.i);
			}
		}
	}
}

class StsObjectFieldTest extends StsObject
{
	int i = 10;
	StsObjectTest parent;

	StsObjectFieldTest()
	{
	}

	StsObjectFieldTest(boolean persistent)
	{
		super(persistent);
	}

	public boolean initialize(StsModel model)
	{
		return true;
	}

	void setParent(StsObjectTest parent)
	{
		this.parent = parent;
	}
}