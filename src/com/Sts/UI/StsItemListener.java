//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.UI;

// This example is from _Java Examples in a Nutshell_. (http://www.oreilly.com)
// Copyright (c) 1997 by David Flanagan
// This example is provided WITHOUT ANY WARRANTY either expressed or implied.
// You may study, use, modify, and distribute it for non-commercial purposes.
// For any commercial use, see http://www.davidflanagan.com/javaexamples

// Modified from StsActionListener to handle an ItemListener. 18 July 98

import com.Sts.Utilities.*;

import java.awt.event.*;
import java.lang.reflect.*;

public class StsItemListener implements ItemListener
{
	protected Object target;
	protected Method m;
	protected Object[] arguments = null;

	/** constructor for an StsItemListener which is a target.method with ItemEvent,String arguments
	 *  @name passed as parameter to invoked method
	 * @param target Object instance object
	 * @methodName String method name
	 */
	public StsItemListener(Object target, String methodName)
	{
		Class c, parameters[];

		this.target = target;

		// get class
		c = target.getClass();

		// get argument classes
		try
		{
			parameters = new Class[] { Class.forName("java.awt.event.ItemEvent") };
		}
		catch(Exception e)
		{
			StsException.outputException("Couldn't find class java.awt.event.ItemEvent", e, StsException.FATAL);
			return;
		}

		// get method
		try
		{
			m = c.getMethod(methodName, parameters);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "constructor", "Couldn't find method: " + methodName + " in class: " + StsToolkit.getSimpleClassname(c), e);
		}
	}

	public StsItemListener(Object target, String methodName, Object[] arguments)
	{
		Class c, parameters[];

		this.target = target;
        this.arguments = arguments;
		// get class
		c = target.getClass();

		// get argument classes
		try
		{
			if(arguments == null)
				parameters = new Class[] { Class.forName("java.awt.event.ItemEvent")};
			else
			{
				parameters = new Class[arguments.length + 1];
				parameters[0] = Class.forName("java.awt.event.ItemEvent");
				for(int n = 0; n < arguments.length; n++)
					parameters[n+1] = arguments[n].getClass();
			}
		}
		catch(Exception e)
		{
			StsException.outputException("Couldn't find class java.awt.event.ItemEvent", e, StsException.FATAL);
			return;
		}

		// get method
		try
		{
			m = c.getMethod(methodName, parameters);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "constructor", "Couldn't find method: " + methodName + " in class: " + StsToolkit.getSimpleClassname(c), e);
		}
	}
	/** constructor for a static method
	 * @className String name of target class
	 * @methodName String name of target method
	 */
	public StsItemListener(String className, String methodName)
	{
		Class c, parameters[];

		this.target = null;

		// get class
		try
		{
			c = Class.forName(className);
		}
		catch(Exception e)
		{
			StsException.outputException("Couldn't find class: " + className,
										 e, StsException.FATAL);
			return;
		}

		// get argument classes
		try
		{
			parameters = new Class[]
				{Class.forName("java.awt.event.ItemEvent")};
		}
		catch(Exception e)
		{
			StsException.outputException("Couldn't find class java.awt.event.ItemEvent", e, StsException.FATAL);
			return;
		}

		// get method
		try
		{
			m = c.getMethod(methodName, parameters);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "constructor", "Couldn't find method: " + methodName + " in class: " + StsToolkit.getSimpleClassname(c), e);
		}
	}

	/** constructor for a static method
	 * @className String name of target class
	 * @methodName String name of target method
	 * @arg Object argument Class(es) to be passed
	 */
	public StsItemListener(Class c, String methodName)
	{
		Class parameters[];

		this.target = null;

		// get argument classes
		try
		{
			parameters = new Class[] {Class.forName("java.awt.event.ItemEvent")};
		}
		catch(Exception e)
		{
			StsException.outputException("Couldn't find class java.awt.event.ItemEvent", e, StsException.FATAL);
			return;
		}

		// get method
		try
		{
			m = c.getMethod(methodName, parameters);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "constructor", "Couldn't find method: " + methodName + " in class: " + StsToolkit.getSimpleClassname(c), e);
		}
	}
 /* TODO: replace target.method with target.getFieldName and target.setFieldName ala beans then classInitialize value with getFieldName so we can support persistent menuItem values

    public static StsItemListener fieldItemListener(Class c, String fieldName)
    {
		if (fieldName == null || c == null) return null;

		try
		{
			get = StsToolkit.getAccessor(c, fieldName, "get", null);
			if(get != null)
			{
				if(!get.isAccessible()) get = null;
			}
			if (get == null)
			{
				try
				{
					field = c.getDeclaredField(fieldName);
					field.setAccessible(true);
				}
				catch (Exception e)
				{
					System.out.println("The field " + fieldName + " in class " + c.getName() + " has no get accessor or doesn't exist.\n");
					editable = false;
				}
			}
			Class returnType;
			if(get != null)
				returnType = get.getReturnType();
			else
				returnType = field.getType();

			set = StsToolkit.getAccessor(c, fieldName, "set", returnType);
			if (set == null && editable)
			{
				System.out.println("The field " + fieldName + " in class " + c.getName() + " has no set accessor or doesn't exist.\n");
				editable = false;
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
*/
    /** When this StsItemListener has its associated itemStateChanged, it calls the target.method with
	 *  arguments of event and name if name is not null, otherwise just the argument event.
	 */

	public void itemStateChanged(ItemEvent event)
	{
		try
		{
			if(arguments != null)
			{
				Object[] allArguments = new Object[arguments.length + 1];
				allArguments[0] = event;
				for(int n = 0; n < arguments.length; n++)
					allArguments[n+1] = arguments[n];
				m.invoke(target, allArguments );
			}
			else
				m.invoke(target, new Object[] {event});
		}
		catch(Exception e)
		{
			StsException.outputException("StsItemListener(itemStateChanged). Couldn't invoke method: "
										 + m.getName() + " in class: " + target.getClass().getName(), e, StsException.FATAL);
		}
	}
}
