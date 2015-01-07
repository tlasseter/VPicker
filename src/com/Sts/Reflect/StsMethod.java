//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Reflect;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.lang.reflect.*;

/**
 * Constructs instance which contains method information
 * for subsequent invocation.
 * Method being invoked must be public in order to find it even if in same package!
 * @author TJLasseter
 */
public class StsMethod
{
    protected Object instance = null;
    protected Object[] methodArgs;
    protected Class c;
    protected Method method;

    /** constructor for an instance method
     * @methodName String method name
     * @methodArgs Object argument Classe(s) to be passed
     */

	public StsMethod()
	{
	}


	public StsMethod(Object instance, String methodName)
	{
		this.instance = instance;
		c = instance.getClass();
		this.methodArgs = methodArgs;
		method = getInstanceMethod(c, methodName, null);
    }

	public StsMethod(Object instance, String methodName, Object methodArg)
    {
        this(instance, methodName, new Object[] {methodArg});
    }

    public StsMethod(Object instance, String methodName, Object[] methodArgs)
    {
        this.instance = instance;
        c = instance.getClass();
        this.methodArgs = methodArgs;
        method = getInstanceMethod(c, methodName, methodArgs);
    }

    public StsMethod(Object instance, String methodName, Class[] argClasses)
    {
        this.instance = instance;
        c = instance.getClass();
 //       this.methodArgs = methodArgs;

        method = getInstanceMethod(c, methodName, argClasses);
    }

	static public void constructInvoke(Object instance, String methodName, Object[] methodArgs)
		throws IllegalAccessException, InvocationTargetException, NoSuchMethodException
	{
		StsMethod method = new StsMethod(instance, methodName, methodArgs);
		method.method.invoke(instance, methodArgs);
	}

    /** Get an instance method
     * @param c Class of instance
     * @methodName String method name
     * @methodArgs Object argument Classe(s) to be passed
     */
    static public Method getInstanceMethod(Class c, String methodName)
    {
        return getInstanceMethod(c, methodName, null);
    }

    static public Method getInstanceMethod(Class c, String methodName, Object[] methodArgs)
    {
        Method method = null;
        Class argClasses[];

        argClasses = getArgClasses(methodArgs);

        // get method
        // Note the use of getMethod here for instance methods and the use of getDeclaredMethod below for static methods

        try
        {
            method = c.getMethod(methodName, argClasses);
            method.setAccessible(true);
        }
        catch (Exception e)
        {
            String argNames = getArgNames(argClasses);
            StsException.outputWarningException(StsMethod.class, "getInstanceMethod(c, methodName, methodArgs)", "Method: " + StsToolkit.getSimpleClassname(c) + "." + methodName + "(" + argNames + ")", e);
        }

        return method;
    }

    static public Method getInstanceMethod(Class c, String methodName, Class[] argClasses)
    {
        Method method = null;

        // get method
        // Note the use of getMethod here for instance methods and the use of getDeclaredMethod below for static methods

        try
        {
			method = c.getMethod(methodName, argClasses);
        }
        catch (Exception e)
        {
            String argNames = getArgNames(argClasses);
            StsException.outputException(errorMessage("getInstanceMethod(c, methodName, argClasses) failed.", c, method, methodName, argClasses),
                                         e, StsException.WARNING);
        }

        return method;
    }
	static private String errorMessage(String message, Class c, Method method, String methodName, Class[] argClasses)
	{
		StringBuffer stringBuffer = new StringBuffer("StsMethod error: " + message);
		if(method != null)
			stringBuffer.append(" methodName: " + method.getName() + " in class " + method.getDeclaringClass().getName());
		else if(c != null)
			stringBuffer.append(" method is null for class " + StsToolkit.getSimpleClassname(c) + " methodName: " + methodName);
		else
			stringBuffer.append(" unknown method in undefined class");
		stringBuffer.append(" method args ");

        if(argClasses != null)
        {
            for(int n = 0; n < argClasses.length; n++)
                stringBuffer.append(" " + argClasses[n].getName());
        }
		else
		    stringBuffer.append(" no class args");

		return stringBuffer.toString();
	}

    static private String errorMessage(String message, Class c, Method method, Object[] methodArgs)
	{
		StringBuffer stringBuffer = new StringBuffer("StsMethod error: " + message);
		if(method != null)
			stringBuffer.append(" methodName: " + method.getName() + " in class " + method.getDeclaringClass().getName());
		else if(c != null)
			stringBuffer.append(" method is null for class " + StsToolkit.getSimpleClassname(c));
		else
			stringBuffer.append(" unknown method in undefined class");
		stringBuffer.append(" method args ");
		if(method != null)
		{
			Class[] argClasses = method.getParameterTypes();
			for(int n = 0; n < argClasses.length; n++)
				stringBuffer.append(" " + argClasses[n].getName());

		}
		else  if(methodArgs != null)
		{
			stringBuffer.append(" method args ");
			for(int n = 0; n < methodArgs.length; n++)
            {
                if(methodArgs[n] != null)
                    stringBuffer.append(" " + methodArgs[n].getClass().getName());
                else
                    stringBuffer.append(" null");
            }
        }
		else
			stringBuffer.append(" none ");

		return stringBuffer.toString();
	}

    public Class getMethodClass()
    {
		return c;
    }

    static private String getArgNames(Class argClasses[])
    {
        if (argClasses == null)
        {
            return new String("None.");
        }

        int length = argClasses.length;
        if (length == 0)
        {
            return new String("None.");
        }

        String string = new String();

        for (int n = 0; n < length; n++)
        {
            String classString = new String(argClasses[n].getName());
            string = new String(string + classString + "\n");
        }

        return string;
    }

    /** Given a list of arguments, get corresponding list of classes
     * so we can identify method with these arguments
     */
    static private Class[] getArgClasses(Object[] methodArgs)
    {
        Class[] argClasses;

        if (methodArgs == null || methodArgs.length == 0 || methodArgs[0] == null)
        {
            argClasses = new Class[0];
        }
        else
        {
            argClasses = new Class[methodArgs.length];
            for (int n = 0; n < methodArgs.length; n++)
            {
                try
                {argClasses[n] = methodArgs[n].getClass();
                }
                catch (Exception e)
                {
                    StsException.outputException(
                        "StsMethod.getArgClasses(Object[]): Couldn't get class for argument Object",
                        e, StsException.FATAL);
                    return null;
                }
            }
        }

        return argClasses;
    }

    /** constructor for a static method from className
     * @className String name of instance class
     * @methodName String name of instance method
     * @methodArgs Object argument Class(es) to be passed
     */

	public StsMethod(String className, String methodName)
	{
		this(className, methodName, new Object[0] );
    }

	public StsMethod(String className, String methodName, Object methodArg)
    {
        this(className, methodName, new Object[] {methodArg});
    }

    public StsMethod(String className, String methodName, Object[] methodArgs)
    {
        try
        {
            c = Class.forName(className);
            this.methodArgs = methodArgs;
            method = getStaticMethod(c, methodName, methodArgs);
        }
        catch (Exception e)
        {
            StsException.outputException(errorMessage("StsMethod constructor(className, methodName, methodArgs) failed.", c, method, methodArgs),
                                         e, StsException.WARNING);
            return;
        }
    }

    public StsMethod(Class c, String methodName, Object methodArg, Class argClass)
    {
        this(c, methodName, new Object[] {methodArg}, new Class[] {argClass} );
    }

    /**
     * Get a static method for Class. Use this method if any of the args are interfaces as reflection cannot find
     * signature of method unless interface class is explicitly defined.
     *
     * @param c class
     * @param methodName name of method
     * @param argClasses list of classes of method arguments
     * @param methodArgs method arguments
     */
    public StsMethod(Class c, String methodName, Object[] methodArgs, Class[] argClasses)
    {
        try
        {
            this.c = c;
            this.methodArgs = methodArgs;
            method = c.getDeclaredMethod(methodName, argClasses);
        }
        catch (Exception e)
        {
            String argNames = getArgNames(argClasses);
            StsException.outputException(errorMessage("StsMethod constructor(methodName, methodArgs, argClasses) failed.", c, method, methodArgs),
                                         e, StsException.WARNING);
        }
    }

    /** Get a static method from Class
     * @param c Class of instance
     * @methodName String method name
     * @methodArgs Object argument Classe(s) to be passed
     */
    static private Method getStaticMethod(Class c, String methodName, Object[] methodArgs)
    {
        Method method = null;
        Class argClasses[];

        argClasses = getArgClasses(methodArgs);

        // get method
        // Note the use of getDeclaredMethod here for static methods and the use of getMethod above for instance methods

        try
        {
            method = c.getDeclaredMethod(methodName, argClasses);
            return method;
        }
        catch (Exception e)
        {
            String argNames = getArgNames(argClasses);

            StsException.outputException(errorMessage("StsMethod.getStaticMethod(c, methodName, methodArgs) failed.", c, method, methodArgs),
                                         e, StsException.WARNING);
            return null;
        }
    }

    /** Get a static method from Class
     * @param c Class of instance
     * @methodName String method name
     * @argClasses Object argument Classe(s) to be passed
     */
    public StsMethod(Class c, String methodName, Class[] argClasses)
    {
        // Note the use of getDeclaredMethod here for static methods and the use of getMethod above for instance methods
        try
        {
            this.c = c;
            method = c.getDeclaredMethod(methodName, argClasses);
        }
        catch (Exception e)
        {
            String argNames = getArgNames(argClasses);

            StsException.outputException(errorMessage("StsMethod.constructor(c, methodName, argClasses) failed.", c, method, null),
										 e, StsException.WARNING);
        }
    }

    /** constructor for a static method from Class
     * @className String name of instance class
     * @methodName String name of instance method
     * @methodArgs Object argument Class(es) to be passed
     */

	public StsMethod(Class c, String methodName)
	{
		this(c, methodName, null);
    }

    public StsMethod(Class c, String methodName, Object methodArg)
    {
        this(c, methodName, new Object[] {methodArg});
    }

    public StsMethod(Class c, String methodName, Object[] methodArgs)
    {
        // get method
        method = getStaticMethod(c, methodName, methodArgs);
        this.methodArgs = methodArgs;
    }

    static public Object invokeStaticMethod(Class c, String methodName, Object[] methodArgs)
    {
		Method method = null;
        try
        {
            StsMethod staticMethod = new StsMethod(c, methodName, methodArgs);
            method = staticMethod.getMethod();
            return method.invoke(null, methodArgs);
        }
        catch (Exception e)
        {
		    StsException.outputException(errorMessage("invokeStaticMethod() failed.\n", c, method, methodArgs), e, StsException.WARNING);
            return null;
        }
    }

    public void invokeInstanceMethod()
    {
        try
        {
            method.invoke(instance, methodArgs);
        }
        catch (Exception e)
        {
            String methodName;

            if (method == null)
            {
                methodName = "null";
            }
            else
            {
                methodName = method.getName();

            }
            StsException.outputException(errorMessage("StsMethod.invokeInstanceMethod() failed.", c, method, methodArgs),
                                         e, StsException.WARNING);
        }
    }

    public void invokeInstanceMethod(Object instance, Object[] methodArgs)
    {
        try
        {
            method.invoke(instance, methodArgs);
        }
        catch (Exception e)
        {
            String methodName;

            if (method == null)
            {
                methodName = "null";
            }
            else
            {
                methodName = method.getName();

            }
            StsException.outputException(errorMessage("StsMethod.invokeInstanceMethod() failed.", c, method, methodArgs),
                                         e, StsException.WARNING);
        }
    }

    public void invokeInstanceMethod(Object[] methodArgs)
    {
        try
        {
            method.invoke(instance, methodArgs);
        }
        catch (Exception e)
        {
            String methodName;

            if (method == null)
            {
                methodName = "null";
            }
            else
            {
                methodName = method.getName();

            }
            StsException.outputException(errorMessage("StsMethod.invokeInstanceMethod() failed.", c, method, methodArgs),
										 e, StsException.WARNING);
        }
    }

    public void invokeStaticMethod(Object[] methodArgs)
    {
        try
        {
            method.invoke(null, methodArgs);
        }
        catch (Exception e)
        {
            StsException.outputException(errorMessage("StsMethod.invokeStaticMethod() failed.", c, method, methodArgs), e, StsException.WARNING);
        }
    }

	/** For beans, we pass in the class and instance methodName; later the instance will be set
	 *  when the beanObject is defined.  So use this constructor to build the instance method.
	 */
	public static StsMethod constructInstanceMethodFromClass(Class c, String methodName)
	{
		StsMethod method = new StsMethod();
		Method m = StsMethod.getInstanceMethod(c, methodName);
		method.method = m;
		return method;
	}


    // Accessors
	public void setInstance(Object instance) { this.instance = instance; }
    public Object getInstance() { return instance; }

    public Object[] getMethodArgs()
    {
		return methodArgs;
    }

    public Method getMethod()
    {
		return method;
    }

	public String toString()
	{
		return this.errorMessage("method ok: ", c, method, methodArgs);
	}

	static public void main(String[] args)
	{
		StsMethod method;
//		StsToolbar toolbar = new StsToolbar("toolbar", false);
//		method = new StsMethod(StsWin3d.class, "closeToolbar", new Class[] { StsToolbar.class } );
//		System.out.println(method.toString());
		StsPreStackVelocityModel velocityModel = new com.Sts.DBTypes.StsPreStackVelocityModel3d();
		method = new StsMethod(velocityModel, "computeInterpolatedVelocityProfile", new Class[] { Integer.class, Integer.class });
		System.out.println(method.toString());

//		StsVelocityAnalysisToolbar velToolbar = new StsVelocityAnalysisToolbar();
//		method = new StsMethod(StsWin3d.class, "closeToolbar", (StsToolbar)velToolbar);
//		System.out.println(method.toString());
	}
}
