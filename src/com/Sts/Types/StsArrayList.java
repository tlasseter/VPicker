package com.Sts.Types;

import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import java.lang.reflect.*;
import java.util.*;

public class StsArrayList extends ArrayList
{

    public StsArrayList()
    {
    }

 	public void forEach(String methodName, Object[] args)
    {
		Object[] list = this.toArray();
        if(list == null) return;

        int nElements = list.length;
        for(int n = 0; n < nElements; n++)
        {
            Object object = list[n];

            try
            {
                StsMethod method = new StsMethod(object, methodName, args);
                Method m = method.getMethod();
        		m.invoke(object, args);
            }
            catch(Exception e)
            {
        	    StsException.outputException("Method: " + methodName + " failed " +
                                             "for object of class: " +
                                             object.getClass().getName(),
                                             e, StsException.WARNING);
            }
        }
    }

 	public void forEach(String methodName)
    {
    	forEach(methodName, null);
    }

    /** If a single argument, make an array of length 1 for compatibility */
 	public void forEach(String methodName, Object arg)
    {
    	forEach(methodName, new Object[]{arg});
    }

    public void addNoRepeat(Object obj)
    {
        if(contains(obj)) return;
        add(obj);
    }

    public void addNoRepeat(Object[] objects)
    {
        if(objects == null) return;
        for(int n = 0; n < objects.length; n++)
            addNoRepeat(objects[n]);
    }
}
