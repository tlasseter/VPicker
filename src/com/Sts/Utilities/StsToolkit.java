package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.UI.*;
import com.Sts.UI.Progress.StsProgressBarDialog;
import com.Sts.UI.Progress.StsProgressPanel;

import javax.jnlp.*;
import javax.swing.*;
import java.awt.*;
import java.beans.*;
import java.io.*;
import java.lang.reflect.*;
import java.net.*;
import java.nio.*;
import java.security.*;
import java.util.*;
import java.util.jar.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version 1.1
 */

public class StsToolkit
{
    static public String[] browsers = new String[]{"iexplore.exe", "firefox.exe", "netscape.exe"};
    static public String[] s2sArchiveTypes = new String[] {"grid.bin", "grid-seismic.bin","LiveWell.bin","well-dev.bin", "well-logs.bin", "well-td.bin", "well-timeLogs", "db.", "Sensor.bin", "sensor.bin","seis","seis3d","seis2d","seisVsp"};
    static final float nullValue = StsParameters.nullValue;

    static public void centerThisOn(Component component, Component parentComponent)
    {
        int parentWidth = parentComponent.getWidth();
        int parentHeight = parentComponent.getHeight();
        Point corner = parentComponent.getLocationOnScreen();
        int centerX = corner.x + parentWidth / 2;
        int centerY = corner.y + parentHeight / 2;
        int width = component.getWidth();
        int height = component.getHeight();
        component.setLocation(centerX - width / 2, centerY - height / 2);
    }

    static public JDialog createDialog(Component component)
    {
        return createDialog(null, component, true, 0, 0);
    }

    static public JDialog createDialog(Component component, boolean modal)
    {
        return createDialog(null, component, modal, 0, 0);
    }

    static public JDialog createDialog(String title, Component component, boolean modal)
    {
        return createDialog(title, component, modal, 0, 0);
    }

    static public JDialog createDialog(Component component, boolean modal, int width, int height)
    {
        return createDialog(null, component, modal, width, height);
    }

    static public JDialog createDialog(String title, Component component, boolean modal, int width, int height)
    {
        JDialog d = new JDialog();
        if(title != null) d.setTitle(title);
        if(width != 0 && height != 0) d.setSize(width, height);
        d.getContentPane().add(component);
        d.setModal(modal);
        d.pack();
        d.setVisible(true);
        return d;
    }

    static public JFrame createJFrame(Component component)
    {
        JFrame frame = new JFrame();
        frame.getContentPane().add(component);
        frame.pack();
        frame.setVisible(true);
        return frame;
    }

    static public JFrame createJFrame(String title, Component component, int width, int height)
    {
        JFrame frame = new JFrame();
        if(title != null) frame.setTitle(title);
        if(width != 0 && height != 0) frame.setSize(width, height);
        frame.getContentPane().add(component);
        frame.pack();
        frame.setVisible(true);
        return frame;
    }

    static public void centerComponentOnFrame(Component component, Component frame)
    {
        int x0 = frame.getX();
        int w0 = frame.getWidth();
        int w1 = component.getWidth();
        int y0 = frame.getY();
        int h0 = frame.getHeight();
        int h1 = component.getHeight();
        int x = x0 + ((w0 - w1) / 2);
        int y = y0 + ((h0 - h1) / 2);
        component.setLocation(x, y);
    }

    static public void centerComponentOnScreen(Component component)
    {
        int width = component.getWidth();
        int height = component.getHeight();
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int screenWidth = screenSize.width;
        int screenHeight = screenSize.height;
        // handle dual screens
        if(screenWidth > 1920 && (width <= screenWidth / 2))
            screenWidth /= 2;
        if(width > screenWidth)
            width = screenWidth;
        if(height > screenHeight)
            height = screenHeight;

        int x = (screenWidth - width) / 2;
        int y = (screenHeight - height) / 2;
        component.setLocation(x, y);
    }

    static public Object deepCopy(Object oldObject)
    {
        ObjectOutputStream oos = null;
        ObjectInputStream ois = null;
        try
        {
            if(oldObject == null)
                return null;
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            oos = new ObjectOutputStream(bos);
            // serialize and pass the object
            oos.writeObject(oldObject);
            oos.flush();
            ByteArrayInputStream bin = new ByteArrayInputStream(bos.toByteArray());
            ois = new ObjectInputStream(bin);
            // return the new object
            return ois.readObject();
        }
        catch(Exception e)
        {
            StsException.outputException("StsDeepCopy.deepCopy() failed.", e, StsException.WARNING);
            return null;
        }
        finally
        {
            try
            {
                oos.close();
                ois.close();
            }
            catch(Exception e)
            {
            }
        }
    }

    /**
     * Copy non-final, non-static, non-transient fields from oldObject to newObject which are of different classes.
     * Do not include fields from superclasses for either.
     */
    static public boolean copyDifferentClasses(Object oldObject, Object newObject)
    {
        Field[] oldFields = getFields(oldObject.getClass());
        TreeMap newFieldsMap = getFieldsTreeMap(newObject.getClass());
        return copyDifferentClasses(oldObject, newObject, oldFields, newFieldsMap);

    }

    /**
     * Copy all non-final and non-static fields from oldObject to newObject which are of different classes.
     * Include fields from superclasses for both.
     * Set the index of the new object to -1 indicating it is nonpersistent and uninitialized.
     */
    static public boolean copyDifferentStsObjectAllFields(StsObject oldObject, StsObject newObject)
    {
        if(!copyDifferentClassesAllFields(oldObject, newObject)) return false;
        newObject.setIndex(-1);
        return true;
    }

    /**
     * Copy all non-final and non-static fields from oldObject to newObject which are of different classes.
     * Include fields from superclasses for both.
     */
    static public boolean copyDifferentClassesAllFields(Object oldObject, Object newObject)
    {
        Field[] oldFields = getAllFields(oldObject.getClass());
        TreeMap newFieldsMap = getAllFieldsTreeMap(newObject.getClass());
        return copyDifferentClasses(oldObject, newObject, oldFields, newFieldsMap);
    }

    /** Copy fields from oldObject to newObject which may be of different classes. */
    static private boolean copyDifferentClasses(Object oldObject, Object newObject, Field[] oldFields, TreeMap newFieldsMap)
    {
        try
        {
            for(int n = 0; n < oldFields.length; n++)
            {
                Class oldType = oldFields[n].getType();
                String oldName = oldFields[n].getName();
                Object oldFieldObject = oldFields[n].get(oldObject);
                if(oldFieldObject == null) continue;
                Field newField = (Field) newFieldsMap.get(oldName);
                if(newField == null) continue;
                if(oldType != newField.getType()) continue;
                newField.set(newObject, oldFieldObject);
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.copyDifferentClasses(oldObject, newObject, oldFields, newFieldsMap) failed for field: ",
                    e, StsException.WARNING);
            return false;
        }
    }

    /** copy fields in subclass thru superclass from oldObject to newObject */
    static public boolean copyToSuperclass(Object oldObject, Object newObject, Class superClass, boolean excludeStsObjects)
    {
        Class subClass = oldObject.getClass();
        return copySubToSuperclass(oldObject, newObject, subClass, superClass, excludeStsObjects);
    }

    /** copy fields in subclass thru superclass from oldObject to newObject */
    static public boolean copySubToSuperclass(Object oldObject, Object newObject, Class superClass)
    {
        if(oldObject == null) return false;
        return copySubToSuperclass(oldObject, newObject, oldObject.getClass(), superClass, true);
    }

    /** copy fields in in class from oldObject to newObject */
    static public boolean copyClass(Object oldObject, Object newObject, Class c)
    {
        if(oldObject == null) return false;
        return copySubToSuperclass(oldObject, newObject, c, c, true);
    }

    /** copy fields in subclass thru superclass from oldObject to newObject */
    static public boolean copySubToSuperclass(Object oldObject, Object newObject, Class subClass, Class superClass, boolean excludeStsObjects)
    {
        Field field = null;
        int modifiers = 0;
        try
        {
            if(oldObject == null || newObject == null) return false;
            TreeMap sortedFields = getFieldsTreeMap(oldObject, subClass, superClass, excludeStsObjects);
            if(sortedFields == null) return false;
            Object[] fields = sortedFields.values().toArray();
            if(fields == null) return false;
            for(int i = 0; i < fields.length; i++)
            {
                field = (Field) fields[i];
                modifiers = field.getModifiers();
                if(Modifier.isFinal(modifiers))
                    continue;
                if(Modifier.isStatic(modifiers))
                    continue;
                field.setAccessible(true);
                Object fieldObject = field.get(oldObject);
                field.set(newObject, fieldObject);
            }
            return true;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " +
                    fieldname + " with modifiers: " + modifiers, e, StsException.WARNING);
            return false;
        }
    }

    static public Object copy(Object oldObject)
    {
        if(oldObject == null) return null;
        try
        {
            Class c = oldObject.getClass();
            Constructor constructor = c.getConstructor(new Class[0]);
            Object copy = constructor.newInstance(new Class[0]);
            if(!copySameClass(oldObject, copy)) return null;
            return copy;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsToolkit.class, "copy(oldObject)", e);
            return null;
        }
    }

    /**
     * Copy fields oldObject to newObject where both are same class.
     * Don't include fields from superclasses.
     */
    static public boolean copy(Object oldObject, Object newObject)
    {
        if(oldObject.getClass() == newObject.getClass())
            return copySameClass(oldObject, newObject);
        else
            return copyDifferentClasses(oldObject, newObject);
    }

    /**
     * Copy fields from oldObject to newObject which are of same class.
     * Don't include fields from superclasses.
     */
    static public boolean copySameClass(Object oldObject, Object newObject)
    {
        Field field = null;
        try
        {
            Field[] fields = getFields(oldObject.getClass());
            for(int i = 0; i < fields.length; i++)
            {
                field = fields[i];
                Object fieldObject = field.get(oldObject);
                field.set(newObject, fieldObject);
            }
            return true;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " + fieldname, e, StsException.WARNING);
            return false;
        }
    }

    static public Object copyObjectFields(Object oldObject)
    {
        return copyObjectFields(oldObject, false);
    }

    /**
     * Return a copy of the non-final, non-static, and non-transient (if includeTransients == false) of this object.
     * Do not include fields from superclasses.
     */
    static public Object copyObjectFields(Object oldObject, boolean includeTransients)
    {
        Field field = null;
        try
        {
            Constructor constructor = oldObject.getClass().getDeclaredConstructor(null);
            Object newObject = constructor.newInstance(null);
            if(!copyObjectFields(oldObject, newObject, includeTransients))
                return null;
            else
                return newObject;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " +
                    fieldname, e, StsException.WARNING);
            return null;
        }
    }

    /**
     * Copy non-final, non-static, and non-transient (if includeTransients == false) from oldObject to new object.
     * Do not include fields from superclasses.
     */
    static public boolean copyObjectFields(Object oldObject, Object newObject, boolean includeTransients)
    {
        Field field = null;
        try
        {

            Field[] fields = getFields(oldObject.getClass(), includeTransients);
            for(int i = 0; i < fields.length; i++)
            {
                field = fields[i];
                Object fieldObject = field.get(oldObject);
                field.set(newObject, fieldObject);
            }
            return true;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " +
                    fieldname, e, StsException.WARNING);
            return false;
        }
    }

    static public Object copyAllObjectFields(Object oldObject)
    {
        return copyAllObjectFields(oldObject, false);
    }

    static public Object copyAllObjectFields(Object oldObject, boolean includeTransients)
    {
        Field field = null;
        try
        {
            Constructor constructor = oldObject.getClass().getDeclaredConstructor(null);
            Object newObject = constructor.newInstance(null);
            if(!copyAllObjectFields(oldObject, newObject, includeTransients))
                return null;
            return newObject;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " +
                    fieldname, e, StsException.WARNING);
            return null;
        }
    }

    static public boolean copyAllObjectFields(Object oldObject, Object newObject, boolean includeTransients)
    {
        Field field = null;
        try
        {

            Field[] fields = getAllFields(oldObject.getClass(), includeTransients);
            for(int i = 0; i < fields.length; i++)
            {
                field = fields[i];
                Object fieldObject = field.get(oldObject);
                // if field in newObject doesn't exist, exception is thrown: ignore
                try
                {
                    field.set(newObject, fieldObject);
                }
                catch(Exception e)
                {
                }
            }
            return true;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " +
                    fieldname, e, StsException.WARNING);
            return false;
        }
    }

    static public boolean copyAllObjectFields(Object oldObject, Object newObject, Class c, boolean includeTransients)
    {
        Field field = null;
        try
        {
            Field[] fields = getAllFields(c, includeTransients);
            for(int i = 0; i < fields.length; i++)
            {
                field = fields[i];
                Object fieldObject = field.get(oldObject);
                field.set(newObject, fieldObject);
            }
            return true;
        }
        catch(Exception e)
        {
            String fieldname = "null";
            if(field != null)
                fieldname = field.getName();
            StsException.outputException("StsToolkit.copy(oldObject, newObject) failed for field: " +
                    fieldname, e, StsException.WARNING);
            return false;
        }
    }

    /**
     * Copy non-final, non-static, and non-transient fields from oldObject to new object
     * Including fields from superclasses.
     */
    static public boolean copyAllObjectNonTransientFields(Object oldObject, Object newObject)
    {
        return copyAllObjectFields(oldObject, newObject, false);
    }

    /**
     * Copy non-final, non-static, and non-transient fields from oldObject to new object.
     * Do not include fields from superclasses.
     */
    static public boolean copyObjectNonTransientFields(Object oldObject, Object newObject)
    {
        return copyObjectFields(oldObject, newObject, false);
    }

    static public Object copyObjectNonTransientFields(Object oldObject)
    {
        return copyAllObjectFields(oldObject, false);
    }

    static public void serializeObject(Object object, String directory, String filename)
    {
        FileOutputStream fos;
        String pathname = directory + File.separator + filename;
        try
        {
            fos = new FileOutputStream(pathname);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(object);
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.serializeObject() failed to write file: " + pathname,
                    e, StsException.WARNING);
            return;
        }
    }

    static public Object deserializeObject(String directory, String filename)
    {
        StsFile file = StsFile.constructor(directory, filename);
        return deserializeObject(file);
    }

    static public boolean doOutputJar(String jarName, String rootDirectory, StsProgressBarDialog progressBarDialog)
    {
        try
        {
            Manifest manifest = new Manifest();
            manifest.getMainAttributes().put(Attributes.Name.MANIFEST_VERSION, "1.0");

            File outFile = new File(rootDirectory + File.separator + "archiveFiles" + File.separator + jarName);
            FileOutputStream bo = new FileOutputStream(outFile);
            JarOutputStream jo = new JarOutputStream(bo, manifest);

            File inFile = new File(rootDirectory);
            addToJar(rootDirectory, inFile, jo, progressBarDialog);

            jo.close();
            bo.close();
            return true;
        }
        catch(Exception ex)
        {
            System.out.println("Failed to write jar file: " + jarName);
            return false;
        }
    }


    static public int numFilesInDirectory(File source)
    {
         BufferedInputStream in = null;
         int cnt = 0;
         try
         {
           if (source.isDirectory())
           {
             for (File nestedFile: source.listFiles())
             {
                 if(nestedFile.isFile())
                 {
                     for(int ii=0; ii<s2sArchiveTypes.length; ii++)
                     {
                        if(nestedFile.getName().contains(s2sArchiveTypes[ii]))
                            cnt++;
                     }
                 }
                 else
                     cnt = cnt + numFilesInDirectory(nestedFile);
             }
             return cnt;
           }
         }
         catch(Exception ex)
         {
            return 0;
         }
         finally
         {
            if (in != null)
               try{in.close();} catch(Exception ex) {}
         }
        return cnt;
    }
    
    static private boolean addToJar(String rootDirectory, File source, JarOutputStream target, StsProgressBarDialog progressDialog)
    {
        BufferedInputStream in = null;
         try
         {
           int len = rootDirectory.length();
           String name = source.getPath().replace("\\", "/");
           if (source.isDirectory())
           {
             if (!name.endsWith("/"))
                 name += "/";
             name = name.substring(len);
             if (!name.isEmpty())
             {
               JarEntry entry = new JarEntry(name);
               entry.setTime(source.lastModified());
               target.putNextEntry(entry);
               target.closeEntry();
             }
             for (File nestedFile: source.listFiles())
             {
                 if(nestedFile.isFile())
                 {
                     for(int ii=0; ii<s2sArchiveTypes.length; ii++)
                     {
                        if(nestedFile.getName().contains(s2sArchiveTypes[ii]))
                        {
                            progressDialog.progressPanel.incrementCount();
                            progressDialog.progressPanel.appendLine("Adding " + nestedFile.getName() + " to archive...");
                            addToJar(rootDirectory, nestedFile, target, progressDialog);
                            break;
                        }
                     }
                 }
                 else
                     addToJar(rootDirectory, nestedFile, target, progressDialog);
             }
             return true;
           }

           name = source.getPath().replace("\\", "/");
           name = name.substring(len);

           JarEntry entry = new JarEntry(name);
           entry.setTime(source.lastModified());
           target.putNextEntry(entry);
           in = new BufferedInputStream(new FileInputStream(source));

           byte[] buffer = new byte[1024];
           while (true)
           {
             int count = in.read(buffer);
             if (count == -1)
               break;
             target.write(buffer, 0, count);
           }
           target.closeEntry();
         }
         catch(Exception ex)
         {
            return false;
         }
         finally
         {
            if (in != null)
               try{in.close();} catch(Exception ex) {}
         }
        return true;
    }
    
    static public Object deserializeObject(StsFile file)
    {

        try
        {
            if(!file.exists())
                return null;
            InputStream in = file.getInputStream();
            ObjectInputStream is = new ObjectInputStream(in);
            return is.readObject();
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.deserializeObject for pathname " + file.getURLDirectory() + File.separator + file.getFilename() + " failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    /**
     * Get fields for this class which are contained in class range subClass thru superclass.
     * Exclude finals, statics, and transients. Also exclude StsObject fields if excludeStsObjects == true.
     */
    static public Object[] getFields(Object object, Class subClass, Class superClass, boolean excludeStsObjects)
    {
        TreeMap sortedMap = getFieldsTreeMap(object, subClass, superClass, excludeStsObjects);
        if(sortedMap == null) return new Object[0];
        return sortedMap.values().toArray();
    }

    /**
     * Get fields TreeMap (key-value pairs sorted by fieldName) for just this class, not superclasses.
     * Exclude finals, statics, and transients.
     */
    static public TreeMap getFieldsTreeMap(Class c)
    {
        return getFieldsTreeMap(c, null, false);
    }

    /**
     * Get fields TreeMap (key-value pairs sorted by fieldName) for this class and superclasses.
     * Exclude finals, statics, and transients.
     */
    static public TreeMap getAllFieldsTreeMap(Class c)
    {
        return getFieldsTreeMap(c, Object.class, false);
    }

    static public TreeMap getFieldsTreeMap(Object object, Class subClass, Class superClass, boolean excludeStsObjects)
    {
        if(object == null) return null;
        if(subClass == null) subClass = object.getClass();
        return getFieldsTreeMap(subClass, superClass, excludeStsObjects);
    }

    /** get all non-transient fields from subClass up thru and including superClass in a TreeMap (key-value pairs sorted by fieldName) */
    static public TreeMap getFieldsTreeMap(Class subClass, Class superClass, boolean excludeStsObjects)
    {
        try
        {
            TreeMap sortedFields = new TreeMap();
            Class nextClass = subClass;

            while (nextClass != null)
            {
                Field[] newFields = getFields(nextClass);
                for(int i = 0; i < newFields.length; i++)
                {
                    if(!excludeStsObjects || !StsObject.class.isAssignableFrom(newFields[i].getType()))
                    {
                        String fieldName = newFields[i].getName();
                        Field currentField = (Field) sortedFields.get(fieldName);
                        if(currentField != null)
                        {
                            StsException.systemError("StsToolkit.getFields() found duplicate field name " + fieldName +
                                    " in classes: " + currentField.getDeclaringClass().getName() +
                                    " and " + newFields[i].getDeclaringClass().getName() + "\n" +
                                    "  Will ignore new field found.");
                        }
                        else
                        {
                            sortedFields.put(fieldName, newFields[i]);
                        }
                    }
                }
                if(superClass == null || nextClass == superClass)
                    nextClass = null;
                else
                    nextClass = nextClass.getSuperclass();
            }
            return sortedFields;
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.getFields() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    /** Get fields for just this class, not superclasses.  Exclude finals, statics, and transients. */
    static private Field[] getFields(Class c)
    {
        return getFields(c, false);
    }

    /**
     * Get fields for just this class, not superclasses.  Exclude finals, statics,
     * and transients (if includeTransients is false)
     */
    static private Field[] getFields(Class c, boolean includeTransients)
    {
        int modifiers = 0;
        try
        {
            Field[] fields = c.getDeclaredFields();
            int nGoodFields = 0;
            Field[] goodFields = new Field[fields.length];
            for(int i = 0; i < fields.length; i++)
            {
                Field field = fields[i];
                modifiers = field.getModifiers();
                if(Modifier.isFinal(modifiers))
                    continue;
                if(Modifier.isStatic(modifiers))
                    continue;
                if(!includeTransients && Modifier.isTransient(modifiers))
                    continue;
                field.setAccessible(true);
                goodFields[nGoodFields++] = field;
            }
            goodFields = (Field[]) StsMath.trimArray(goodFields, nGoodFields);
            return goodFields;
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.getFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public Field getField(String classname, String fieldname)
    {
        Class classClass = getClassForName(classname);
        if(classClass == null)
        {
            StsException.systemError(StsToolkit.class, "getField", " can't find class " + classname);
            return null;
        }
        return getField(classClass, fieldname);
    }

    static public Field getField(Class classClass, String fieldname)
    {
        while (classClass != null)
        {
            try
            {
                Field field = classClass.getDeclaredField(fieldname);
                if(field != null)
                {
                    field.setAccessible(true);
                    return field;
                }
            }
            catch(NoSuchFieldException nsfe)
            {
            }
            catch(Exception e)
            {
                StsException.outputWarningException(StsToolkit.class, "getField", e);
                return null;
            }
            classClass = classClass.getSuperclass();
        }
        StsException.systemError(StsToolkit.class, "getField", " can't find field " + fieldname + " in class " + classClass.getName());
        return null;
    }

    static public Class getClassForName(String className)
    {
        try
        {
            return StsToolkit.class.getClassLoader().loadClass(className);
        }
        catch(Exception e)
        {
            return null;
        }
    }

    /** Get all non-final, non-static, and non-transient for this class including superclasses. */
    static public Field[] getAllFields(Class c)
    {
        return getAllFields(c, false);
    }

    /**
     * Get all non-final, non-static, and non-transient (if includeTransients == false) for this class
     * including superclasses.
     */
    static public Field[] getAllFields(Class c, boolean includeTransients)
    {
        Field[] fields = new Field[0];
        while (c != null)
        {
            Field[] newFields = getFields(c, includeTransients);
            fields = (Field[]) StsMath.arrayAddArray(fields, newFields);
            c = c.getSuperclass();
        }
        return fields;
    }

    static public void objectToString(String name, Object object)
    {
        System.out.println(name);
        objectToPrint(System.out, object, null, null, "");
        System.out.println("");
    }

    static public void objectToPrint(PrintStream printStream, String name, Object object, Class superClass)
    {
        System.out.println(name);
        objectToPrint(printStream, object, null, superClass, "");
        System.out.println("");
    }

    static public void objectToPrint(PrintStream printStream, String name, Object object, Class subClass, Class superClass)
    {
        System.out.println(name);
        objectToPrint(printStream, object, subClass, superClass, "");
        System.out.println("");
    }

    static private void objectToPrint(PrintStream printStream, Object object, Class subClass, Class superClass, String prefix)
    {
        try
        {
            if(subClass == null) subClass = object.getClass();

            Object[] fields = StsToolkit.getFields(object, subClass, superClass, false);
            int nFields = fields.length;
            for(int n = 0; n < nFields; n++)
            {
                Field field = (Field) fields[n];
                if(field.get(object) == null) continue;
                if(Modifier.isTransient(field.getModifiers())) continue;

                String fieldName = field.getName();
                Class type = field.getType();
                if(type.isPrimitive() || type == String.class)
                    printStream.print(" " + prefix + fieldName + " " + field.get(object).toString());
                else if(type.isArray())
                {
                    printStream.print(" " + prefix + fieldName);

                    Class arrayType = type.getComponentType();
                    if(arrayType.isArray())
                    {
                        StsException.systemError("Developer!  StsParameterFile cannot handle multidimensioned arrays for field " + fieldName + " in object " + StsToolkit.getSimpleClassname(object) + "!  Sorry!");
                        return;
                    }
                    if(arrayType.equals(Integer.TYPE))
                    {
                        int[] ints = (int[]) field.get(object);
                        for(int i = 0; i < ints.length; i++)
                            printStream.print(" " + ints[i]);
                    }
                    else if(arrayType.equals(Float.TYPE))
                    {
                        float[] floats = (float[]) field.get(object);
                        for(int i = 0; i < floats.length; i++)
                            printStream.print(" " + floats[i]);
                    }
                    else if(arrayType.equals(Double.TYPE))
                    {
                        double[] doubles = (double[]) field.get(object);
                        for(int i = 0; i < doubles.length; i++)
                            printStream.print(" " + doubles[i]);
                    }
                    else if(arrayType.equals(Boolean.TYPE))
                    {
                        boolean[] booleans = (boolean[]) field.get(object);
                        for(int i = 0; i < booleans.length; i++)
                        {
                            if(booleans[i])
                                printStream.print(" true");
                            else
                                printStream.print(" false");
                        }
                    }
                    else if(arrayType.equals(Byte.TYPE))
                    {
                        byte[] bytes = (byte[]) field.get(object);
                        for(int i = 0; i < bytes.length; i++)
                            printStream.print(" " + bytes[i]);
                    }
                    else if(arrayType.equals(Character.TYPE))
                    {
                        char[] chars = (char[]) field.get(object);
                        for(int i = 0; i < chars.length; i++)
                            printStream.print(" " + chars[i]);
                    }
                    else if(arrayType.equals(Short.TYPE))
                    {
                        short[] shorts = (short[]) field.get(object);
                        for(int i = 0; i < shorts.length; i++)
                            printStream.print(" " + shorts[i]);
                    }
                    else if(arrayType.equals(Long.TYPE))
                    {
                        long[] longs = (long[]) field.get(object);
                        for(int i = 0; i < longs.length; i++)
                            printStream.print(" " + longs[i]);
                    }
                    else if(arrayType.equals(String.class))
                    {
                        String[] strings = (String[]) field.get(object);
                        for(int i = 0; i < strings.length; i++)
                        {
                            strings[i].replaceAll(" ", "");
                            printStream.print(" " + strings[i]);
                        }
                    }
                    else
                    {
                        Object fieldObject = field.get(object);
                        if(fieldObject != null && fieldObject instanceof Object[])
                        {
                            Object[] objects = (Object[]) fieldObject;
                            int nObjects = objects.length;
                            printStream.print(nObjects + " ");
                            for(int i = 0; i < nObjects; i++)
                                objectToPrint(printStream, objects[i], type, type, prefix + fieldName + ".");
                        }
                    }
                    printStream.print("\n");
                }
                else
                {
                    objectToPrint(printStream, field.get(object), null, null, prefix + fieldName + ".");
                }
            }
        }
        catch(Exception e)
        {
            printStream.print(" none");
        }
    }

    static public void catchUp()
    {
        System.gc();
        System.runFinalization();
        for(int i = 0; i < 10; i++)
            if(Toolkit.getDefaultToolkit().getSystemEventQueue().peekEvent() != null)
            {
                try
                {
                    Thread.sleep(100);
                }
                catch(Exception e)
                {
                }
            }
    }

    static public void sleep(long millisecs)
    {
        try
        {
            Thread.sleep(millisecs);
        }
        catch(Exception e)
        {
        }
    }

    /** WriteObjectXML serializes an instance of Object as XML. Object must implement Serializable */
    public static void writeObjectXML(Object object, String filename)
    {
        writeObjectXML(object, filename, false);
    }

    public static void writeObjectXML(Object object, String filename, boolean append)
    {
        FileOutputStream fstream = null;
        XMLEncoder ostream = null;

        try
        {
            fstream = new FileOutputStream(filename, append);
            ostream = new XMLEncoder(fstream);
            ostream.writeObject(object);
            ostream.flush();
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.writeObjectXML() failed for object: " +
                    object.toString() + " filename: " + filename + ".",
                    e, StsException.WARNING);
        }
        finally
        {
            try
            {
                if(ostream != null)
                    ostream.close();
                if(fstream != null)
                    fstream.close();
            }
            catch(Exception e)
            {
            }
        }
    }

    /** readObjectXML reads an XML file and creates an object. */
    public static Object readObjectXML(String filename) throws Exception
    {
        FileInputStream fstream = null;
        XMLDecoder istream = null;

        try
        {
            fstream = new FileInputStream(filename);
            istream = new XMLDecoder(fstream);
            Object obj = istream.readObject();
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.writeObjectXML() failed for filename: " + filename + ".",
                    e, StsException.WARNING);
        }
        finally
        {
            try
            {
                if(fstream != null)
                    fstream.close();
                if(istream != null)
                    istream.close();
            }
            catch(Exception e)
            {
            }
            return null;
        }
    }

    /** Get class stemnames from classpath and package name.  They may be in a directory of classes or in a jar file. */
    static public ArrayList findFilePackageClasses(String packageName)
    {

        ArrayList classes = new ArrayList();
        boolean debug = false;

        StringBuffer buffer = new StringBuffer(System.getProperty("java.class.path"));

        try
        {
            if(debug)
                System.out.println("StsToolkit.findPackageClasses(): ");
            StringTokenizer pathTokens = new StringTokenizer(buffer.toString(), File.pathSeparator);

            while (pathTokens.hasMoreElements())
            {
                String parentName = (String) pathTokens.nextElement();
                File parentFile = new File(parentName);
                if(parentFile.isDirectory()) // file is a directory; list files in package name subdirectory
                {
                    String packageDirname = packageName.replace('.', File.separatorChar);
                    if(debug)
                        System.out.println("    looking for class files in directory: " + parentName);
                    File packageDirectory = new File(parentFile, packageDirname);
                    if(packageDirectory.isDirectory())
                    {
                        String[] list = packageDirectory.list();
                        for(int i = 0; i < list.length; ++i)
                        {
                            if(!list[i].endsWith(".class"))
                                continue;
                            String classStemName = list[i].replaceAll(".class", "");
                            classes.add(classStemName);
                        }
                    }
                }
                else // file is a jar file: get classes in package
                {
                    String packageDirname = packageName.replace('.', '/');
                    if(!packageDirname.endsWith("/"))
                        packageDirname = packageDirname + "/";
                    if(debug)
                        System.out.println("    looking for class files in jar: " + parentName);
                    JarFile jarfile = new JarFile(parentFile);
                    Enumeration enumeration = jarfile.entries();
                    while (enumeration.hasMoreElements())
                    {
                        JarEntry jarEntry = (JarEntry) enumeration.nextElement();
                        String entryName = jarEntry.getName();
                        String classStemName = truncateClassName(entryName, packageDirname);
                        if(classStemName != null)
                            classes.add(classStemName);
                    }
                }
            }
            return classes;
        }
        catch(Exception e)
        {
            StsException.systemError("StsToolkit.findPackageClasses() failed for package: " + packageName);
            return classes;
        }
    }

    static public String[] getTokens(String string)
    {
        StringTokenizer tokenizer = new StringTokenizer(string);
        int nTokens = tokenizer.countTokens();
        String[] tokens = new String[nTokens];
        int n = 0;
        while (tokenizer.hasMoreTokens())
            tokens[n++] = tokenizer.nextToken();
        return tokens;
    }

    static public boolean launchWebPageViaBrowser(String page)
    {
        boolean result = false;
        for(int i = 0; i < browsers.length; i++)
        {
            result = launchWebPage(page, browsers[i]);
            if(result)
                return true;
        }
        return false;
    }

    static public boolean launchWebPage(String page)
    {
        DownloadService downloadService = JNLPUtilities.getDownloadService();
        if(downloadService == null)
        {
            return launchWebPageViaBrowser(page);
        }
        if(JNLPUtilities.getBasicService().isWebBrowserSupported())
        {
            try
            {
                URL website = new URL(page);
                JNLPUtilities.getBasicService().showDocument(website);
                return true;
            }
            catch(Exception g)
            {
                return launchWebPageViaBrowser(page);
            }
        }
        return true;
    }

    static public boolean launchWebPage(String page, String browser)
    {
        try
        {
            String executeCmd = browser + " " + page;
            Runtime rt = Runtime.getRuntime();
            rt.exec(executeCmd);
            return true;
        }
        catch(Exception g)
        {
            return false;
        }
    }

    static public String truncateClassName(String className, String packageDirname)
    {
        int start = className.indexOf(packageDirname);
        if(start < 0)
            return null;
        int end = className.indexOf(".class");
        if(end < 0)
            return null;
        start += packageDirname.length();
        return className.substring(start, end);
    }

    /**
     * Given a filename, for example, "com/Sts/Workflow/PlugIn/MyPlugIn.class" and a package name such as
     * "com.Sts.Workflow.PlugIn", return "MyPlugIn" since the package (with '.' replaced by '/') matches.
     *
     * @param filenames   Names of files we wish to check if they match package names. Use "/" as separator.
     * @param packageName Name of package we wish to match.
     * @return matching class names
     */
    static public ArrayList truncateClassnames(String[] filenames, String packageName)
    {
        ArrayList classes = new ArrayList();

        String packageDirname = packageName.replace('.', '/');
        if(!packageDirname.endsWith("/"))
            packageDirname = packageDirname + "/";

        System.out.println("truncateClassNames() packageName: " + packageName);
        for(int n = 0; n < filenames.length; n++)
        {
            String className = filenames[n];
            System.out.println("                     filename: " + className);
            int start = className.indexOf(packageDirname);
            if(start < 0)
                continue;
            int end = className.indexOf(".class");
            if(end < 0)
                continue;
            start += packageDirname.length();
            className = className.substring(start, end);
            className.replace('/', '.');
            System.out.println("                     className: " + className);
            classes.add(className);
        }
        return classes;
    }

    static public String toDebugString(Object object)
    {
        return getSimpleClassname(object) + "@" + Integer.toHexString(object.hashCode());
    }

    static public Method getAccessorX(Class c, String fieldName, String prefix, Class arg)
    {
        String cap = fieldName.substring(0, 1);
        String accessorName = prefix + cap.toUpperCase() + fieldName.substring(1);
        Method accessor = null;
        Class iclass = c;
        for(; iclass != null && accessor == null; iclass = iclass.getSuperclass())
        {
            try
            {
                Class[] argClass = (arg == null ? null : new Class[]{arg});
                accessor = iclass.getDeclaredMethod(accessorName, argClass);
                accessor.setAccessible(true);
            }
            catch(Exception e)
            {
                if(arg != null && arg.equals(Boolean.TYPE))
                {
                    Class[] argClass = new Class[]
                            {Boolean.TYPE};
                    try
                    {
                        accessor = iclass.getDeclaredMethod(accessorName, argClass);
                    }
                    catch(Exception ex)
                    {
                    }
                }
            }
        }
        return accessor;
    }

	static public Method getAccessor(Class c, String fieldName, String prefix, Class arg)
	{
		String cap = fieldName.substring(0, 1);
		String accessorName = prefix + cap.toUpperCase() + fieldName.substring(1);
		Class iclass = c;
		for(; iclass != null; iclass = iclass.getSuperclass())
		{
			try
			{
				Method[] methods = iclass.getDeclaredMethods();
				if(methods == null) continue;
				for(Method accessor : methods)
				{
					if(!accessor.getName().equals(accessorName)) continue;
					accessor.setAccessible(true);
					if(arg == null) return accessor;
					Class[] argTypes = accessor.getParameterTypes();
					if(argTypes.length > 1) continue;
					if(argTypes[0] != arg) continue;
					return accessor;
				}
			}
			catch(Exception e)
			{
				if(arg != null && arg.equals(Boolean.TYPE))
				{
					Class[] argClass = new Class[] {Boolean.TYPE};
					try
					{
						Method accessor = iclass.getDeclaredMethod(accessorName, argClass);
						if(accessor != null) return accessor;
					}
					catch(Exception ex)
					{
						StsException.outputWarningException(StsToolkit.class, "getAccessor", "Failed. couldn't find class.method " + c.getName() + "." + fieldName +
						" with args " + Arrays.toString(argClass), ex);
					}
				}
				return null;
			}
		}
		return null;
	}

    public static boolean clean(final ByteBuffer buffer)
    {
        if(buffer == null || !buffer.isDirect())
        {
            return false;
        }
        Boolean b = (Boolean) AccessController.doPrivileged
                (
                        new PrivilegedAction()
                        {
                            public Object run()
                            {
                                Boolean success = Boolean.FALSE;
                                try
                                {
                                    Method getCleanerMethod = buffer.getClass().getMethod("cleaner", (Class[]) null);
                                    getCleanerMethod.setAccessible(true);
                                    Object cleaner = getCleanerMethod.invoke(buffer, (Object[]) null);
                                    Method clean = cleaner.getClass().getMethod("clean", (Class[]) null);
                                    clean.invoke(cleaner, (Object[]) null);
                                    success = Boolean.TRUE;
                                }
                                catch(Exception e)
                                {
                                    //						StsException.outputException("StsSeismicVolume.clean() failed.",
                                    //													 e, StsException.WARNING);
                                    ;
                                }
                                return success;
                            }
                        }
                );
        return b.booleanValue();
    }

    static public void beep()
    {
        Toolkit.getDefaultToolkit().beep();
    }

    /**
     * Compares two objects.  If they are null or of different classes, returns false.
     * If fields are primitives, and they don't compare return false.
     * If fields are objects, recursively compare objects.
     */

    static public boolean compareObjects(Object object1, Object object2)
    {
        if(object1 == null && object2 == null)
            return true;
        if(object1 == object2)
            return true;
        if(object1 == null || object2 == null)
            return false;
        if(object1.getClass() != object2.getClass())
            return false;
        try
        {
            Field[] fields = getAllFields(object1.getClass(), false);
            for(int n = 0; n < fields.length; n++)
            {
                Field field = fields[n];
                int modifiers = field.getModifiers();
                if(Modifier.isFinal(modifiers))
                    continue;
                if(Modifier.isStatic(modifiers))
                    continue;
                field.setAccessible(true);
                Object fieldObject1 = field.get(object1);
                Object fieldObject2 = field.get(object2);
                Class type = field.getType();
                if(!type.isPrimitive())
                {
                    if(!compareObjects(fieldObject1, fieldObject2))
                        return false;
                }
                else
                {
                    if(!fieldObject1.equals(fieldObject2))
                        return false;
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.compareObjects() failed.", e, StsException.WARNING);
            return false;
        }
    }

    static public boolean compareObjects(Object object1, Object object2, boolean includeTransients)
    {
        if(object1 == null || object2 == null)
            return false;
        if(object1.getClass() != object2.getClass())
            return false;
        try
        {
            Field[] fields = getFields(object1.getClass(), includeTransients);
            for(int n = 0; n < fields.length; n++)
            {
                Field field = fields[n];
                Object fieldObject1 = field.get(object1);
                Object fieldObject2 = field.get(object2);
                Class type = field.getType();
                if(!type.isPrimitive())
                {
                    if(!compareObjects(fieldObject1, fieldObject2))
                        return false;
                }
                else
                {
                    if(!fieldObject1.equals(fieldObject2))
                        return false;
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsToolkit.compareObjects() failed.", e, StsException.WARNING);
            return false;
        }
    }

    static public Icon getCloseIcon()
    {
        return UIManager.getIcon("InternalFrame.closeIcon");
    }

    static public void accumulateHistogram(int[] dataCnt, float dataMin, float dataMax, float value)
    {
        int nSamples = 0;
        accumulateHistogram(dataCnt, dataMin, dataMax, nSamples, value);
    }

    static public void accumulateHistogram(int[] dataCnt, float dataMin, float dataMax, int ttlHistogramSamples, float value)
    {
        float scaledFloat = 254 * (value - dataMin) / (dataMax - dataMin);
        int scaledInt = StsMath.minMax(Math.round(scaledFloat), 0, 254);
        dataCnt[scaledInt] += 1;
        ttlHistogramSamples++;
    }

    static public float[] calculateHistogram(int[] dataCnt, int ttlHistogramSamples)
    {
        float[] dataHist = new float[255];
        for(int i = 0; i < 255; i++)
            dataHist[i] = ((float) dataCnt[i] / (float) ttlHistogramSamples) * 100.0f;
        return dataHist;
    }

    static public float[] buildHistogram(byte[] data, float min, float max)
    {
        int dataCnt[] = new int[255];
        int ttlSamples = data.length;

        if(data == null)
            return null;
        for(int i = 0; i < data.length; i++)
        {
            int indx = StsMath.signedByteToUnsignedInt(data[i]);
            if(indx > 254) continue;
            (dataCnt[indx])++;
        }
        return calculateHistogram(dataCnt, ttlSamples);
    }

    static public float[] buildHistogram(float[] data, float min, float max)
    {
        int dataCnt[] = new int[255];
        int ttlSamples = data.length;

        if(data == null)
            return null;
        for(int i = 0; i < data.length; i++)
            StsToolkit.accumulateHistogram(dataCnt, min, max, data[i]);
        return calculateHistogram(dataCnt, ttlSamples);
    }

    static public float[] buildHistogram(float[][] data, float min, float max)
    {
        int dataCnt[] = new int[255];
        int ttlSamples = data.length;

        if(data == null)
            return null;
        for(int j = 0; j < data.length; j++)
        {
            ttlSamples = +ttlSamples + data[j].length;
            for(int i = 0; i < data[j].length; i++)
                StsToolkit.accumulateHistogram(dataCnt, min, max, data[j][i]);
        }
        return calculateHistogram(dataCnt, ttlSamples);
    }

    static public boolean newDataInDirectory(String source, String[] filterStrings, long lastTime)
    {
        boolean newData = false;
        File dirFile = new File(source);
        if(!dirFile.isAbsolute() || (!dirFile.isDirectory() && !dirFile.isFile()))
            return false;
        StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);
        StsFileSet fileSet = StsFileSet.constructor(source, filter);
        StsAbstractFile[] files = fileSet.getFiles();
        for(int i = 0; i < files.length; i++)
        {
            if(lastTime < files[i].getFile().lastModified())
                return true;
        }
        return false;
    }

    static public long getFileSize(String source)
    {
        File inFile = new File(source);
        if(!inFile.exists()) return 0;
        return inFile.length();
    }

    static public boolean newData(String source, long lastTime)
    {
        File dirFile = new File(source);
        if(!dirFile.isAbsolute() || (!dirFile.isDirectory() && !dirFile.isFile()))
            return false;
        if(lastTime < dirFile.lastModified())
            return true;

        return false;
    }

    /**
     * Create a worker thread for this runnable and run it unless the current thread is already a worker thread
     * (not EventQueue or Java2d.flusher threads).
     *
     * @param runnable new runnable to to be run
     */
    static public void runRunnable(Runnable runnable)
    {
//         String name = Thread.currentThread().getName();
//         System.out.println("runRunnable for thread: " + name);

        if(SwingUtilities.isEventDispatchThread())
        // if (Java2D.isQueueFlusherThread() || SwingUtilities.isEventDispatchThread())
        {
            Thread thread = new Thread(runnable);
            thread.start();
        }
        else
        {
            runnable.run();
        }
    }

    static public Thread checkRunRunnable(Runnable runnable)
    {
//         String name = Thread.currentThread().getName();
//         System.out.println("runRunnable for thread: " + name);

        if(SwingUtilities.isEventDispatchThread())
        // if (Java2D.isQueueFlusherThread() || SwingUtilities.isEventDispatchThread())
        {
            Thread thread = new Thread(runnable);
            thread.start();
            return thread;
        }
        else
        {
            runnable.run();
            return Thread.currentThread();
        }
    }

    /**
     * Run this runnable on the event thread if it is on a worker thread or on the Java2d.flusher thread;
     * if already on eventThread then run it on eventThread.
     *
     * @param runnable
     */
    static public void runLaterOnEventThread(Runnable runnable)
    {
        if(SwingUtilities.isEventDispatchThread())
        {
            // System.out.println("RunLaterOnEventThread. Already on event thread. runnable: " + runnable.toString());
            runnable.run();
        }
        else
        {
            try
            {
                // System.out.println("RunLaterOnEventThread. runnable: " + runnable.toString());
                SwingUtilities.invokeLater(runnable);
            }
            catch(Exception e)
            {
                StsException.outputException("StsToolkit.runLaterOnEventThread(Runnable) failed.", e, StsException.WARNING);
            }
        }
    }

    static public void paintImmediately(Runnable runnable, JComponent component)
    {
        if(SwingUtilities.isEventDispatchThread())
        {
            runnable.run();
            Rectangle r = new Rectangle();
            component.getBounds(r);
            component.paintImmediately(r);
        }
        else
        {
            try
            {
                SwingUtilities.invokeLater(runnable);
            }
            catch(Exception e)
            {
                StsException.outputException("StsToolkit.runLaterOnEventThread(Runnable) failed.", e, StsException.WARNING);
            }
        }
    }

    /**
     * If on a worker thread, invoke and wait on the EventThread.  If already on the EventThread, continue running it.
     * If on the Java2D.flusher thread, run later on the EventThread, though this situation should be flagged as a
     * potential problem.
     *
     * @param runnable
     */
    static public void runWaitOnEventThread(Runnable runnable)
    {
        if(SwingUtilities.isEventDispatchThread())
        {
            // System.out.println("runWaitOnEventThread. Already on event thread runnable: " + runnable.toString());
            runnable.run();
        }
        /*
        else if(Java2D.isOGLPipelineActive() && Java2D.isQueueFlusherThread())
        {
            try
            {
                StsException.systemError(StsToolkit.class, "RunWaitOnEventThread called on Java2D flusher thread: " + runnable.toString());
                SwingUtilities.invokeLater(runnable);
            }
            catch (Exception e)
            {
                StsException.outputException("StsToolkit.runWaitOnEventThread(Runnable) from Java2D flusher thread failed.", e, StsException.WARNING);
            }
        }
        */
        else
        {
            try
            {
                SwingUtilities.invokeAndWait(runnable);
            }
            catch(Exception e)
            {
                StsException.outputException("StsToolkit.runWaitOnEventThread(Runnable) failed.", e, StsException.WARNING);
            }
        }
    }

    public static String getSimpleClassname(Object object)
    {
        if(object == null) return "Null";
        return getSimpleClassname(object.getClass());
    }

    public static String getSimpleClassname(Class c)
    {
        if(c == null) return "Null";
        String classname = c.getName();
        int lastDotIndex = classname.lastIndexOf('.') + 1;
        return classname.substring(lastDotIndex);
    }

    static public String findDirectoryForFile(Frame frame, String directory, String filename)
    {
        try
        {
            String title = "Find directory for file: " + filename;
            StsDirectorySelectionDialog findDirectoryDialog = new StsDirectorySelectionDialog(frame, title, directory, filename, true);
            findDirectoryDialog.setVisible(true);
            directory = findDirectoryDialog.getCurrentDirectory();
            if(directory == null) return null;
            return directory + File.separator;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsToolkit.class, "browseForDirectory", e);
            return null;
        }
    }

    static public final byte[] constructByteArray(int length, byte value)
    {
        byte[] bytes = new byte[length];
        for(int n = 0; n < length; n++)
            bytes[n] = value;
        return bytes;
    }

    static public float[][][] computeSmoothNormals(float[][] pointsZ, int nRows, int nCols, float xInc, float yInc)
    {
        int i, j;

        float[][][] normals = new float[nRows][nCols][];

        float[] idif = new float[3];
        float[] jdif = new float[3];

        idif[1] = yInc;
        jdif[0] = xInc;

        for(i = 0; i < nRows; i++)
        {
            for(j = 0; j < nCols; j++)
            {
                idif[2] = getRowDif(pointsZ, nRows, nCols, i, j);
                jdif[2] = getColDif(pointsZ, nRows, nCols, i, j);
                normals[i][j] = StsMath.crossProduct(idif, jdif);
            }
        }
        return normals;
    }

    static final private float getRowDif(float[][] pointsZ, int nRows, int nCols, int row, int col)
    {
        if(nRows == 1) return 0.0f;
        int rowMin = Math.max(row - 1, 0);
        float zMin = pointsZ[rowMin][col];
        if(zMin == nullValue && rowMin < row)
        {
            rowMin++;
            zMin = pointsZ[rowMin][col];
        }
        if(zMin == nullValue) return 0.0f;
        int rowMax = Math.min(row + 1, nRows - 1);
        float zMax = pointsZ[rowMax][col];
        if(zMax == nullValue && rowMax > row)
        {
            rowMax--;
            zMax = pointsZ[rowMax][col];
        }
        if(zMax == nullValue) return 0.0f;
        int dRow = rowMax - rowMin;
        if(dRow == 0)
            return 0.0f;
        if(dRow == 1)
            return zMax - zMin;
        else
            return (zMax - zMin) / 2;
    }

    static final private float getColDif(float[][] pointsZ, int nRows, int nCols, int row, int col)
    {
        if(nCols == 1) return 0.0f;
        int colMin = Math.max(col - 1, 0);
        float zMin = pointsZ[row][colMin];
        if(zMin == nullValue && colMin < col)
        {
            colMin++;
            zMin = pointsZ[row][colMin];
        }
        if(zMin == nullValue) return 0.0f;
        int colMax = Math.min(col + 1, nCols - 1);
        float zMax = pointsZ[row][colMax];
        if(zMax == nullValue && colMax > col)
        {
            colMax--;
            zMax = pointsZ[row][colMax];
        }
        if(zMax == nullValue) return 0.0f;
        int dRow = colMax - colMin;
        if(dRow == 0)
            return 0.0f;
        if(dRow == 1)
            return zMax - zMin;
        else
            return (zMax - zMin) / 2;
    }

    static public float[][][] computeSmoothNormalsOld(float[][] pointsZ, int nRows, int nCols, float xInc, float yInc)
    {
        int i, j;

        float[][][] normals = new float[nRows][nCols][3];

        float[] idif = new float[3];
        float[] jdif = new float[3];

        idif[1] = yInc;
        jdif[0] = xInc;

        for(i = 0; i < nRows - 1; i++)
        {
            for(j = 0; j < nCols - 1; j++)
            {
                if(pointsZ[i][j] == nullValue) continue;
                if(pointsZ[i + 1][j] == nullValue) continue;
                if(pointsZ[i][j + 1] == nullValue) continue;
                idif[2] = pointsZ[i + 1][j] - pointsZ[i][j];
                jdif[2] = pointsZ[i][j + 1] - pointsZ[i][j];
                float[] xp = StsMath.crossProduct(idif, jdif);
                StsMath.vectorAdd(normals[i][j], xp);
                StsMath.vectorAdd(normals[i + 1][j], xp);
                StsMath.vectorAdd(normals[i][j + 1], xp);
            }
            // Copy normal for last column from second to last column
            normals[i][nCols - 1] = normals[i][nCols - 2];
        }
        if(nRows < 2) return normals;
        // Copy normals for top row from second to top row
        for(j = 0; j < nCols; j++)
        {
            normals[nRows - 1][j] = normals[nRows - 2][j];
        }
        return normals;
    }

    static public float[][][] computeFaceNormals(float[][] pointsZ, int nRows, int nCols, float xInc, float yInc)
    {
        int i, j;

        float[][][] normals = new float[nRows][nCols][3];

        float[] idif = new float[3];
        float[] jdif = new float[3];

        idif[1] = yInc;
        jdif[0] = xInc;

        for(i = 1; i < nRows; i++)
        {
            for(j = 1; j < nCols; j++)
            {
                idif[2] = pointsZ[i][j] - pointsZ[i - 1][j];
                jdif[2] = pointsZ[i][j] - pointsZ[i][j - 1];
                normals[i][j] = StsMath.crossProduct(idif, jdif);
            }
            // Copy normal for 1st column from second
            normals[i][0] = normals[i][1];
        }

        // Copy normals for top row from second to top row
        for(j = 0; j < nCols; j++)
        {
            normals[0][j] = normals[1][j];
        }
        return normals;
    }

    static public float[][][] computeGridCenterNormals(float[][][] normals, int nRows, int nCols, int length, int min)
    {
        float[][][] centerNormals = new float[nRows][nCols][3];

        for(int row = 0; row < nRows - 1; row++)
            for(int col = 0; col < nCols - 1; col++)
                centerNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row + 1][col], normals[row][col + 1], normals[row + 1][col + 1], length, min);
        return centerNormals;
    }

    static public float[][] computeGridCenterValues(float[][] values, int nRows, int nCols, float nullValue)
    {
        float[][] centerValues = new float[nRows][nCols];

        for(int row = 0; row < nRows - 1; row++)
        {
            for(int col = 0; col < nCols - 1; col++)
                centerValues[row][col] = StsMath.average(values[row][col], values[row + 1][col], values[row][col + 1], values[row + 1][col + 1], nullValue);
            centerValues[row][nCols - 1] = nullValue;
        }
        Arrays.fill(centerValues[nRows - 1], nullValue);
        return centerValues;
    }

    static public void print(double[] values, String message)
    {
        System.out.print(message + ": ");
        for(int n = 0; n < values.length; n++)
            System.out.print(" " + values[n]);
        System.out.println();
    }

    static public void print(float[] values, String message)
    {
        System.out.print(message + ": ");
        for(int n = 0; n < values.length; n++)
            System.out.print(" " + values[n]);
        System.out.println();
    }

    static public void print(int[] values, String message)
    {
        System.out.print(message + ": ");
        for(int n = 0; n < values.length; n++)
            System.out.print(" " + values[n]);
        System.out.println();
    }

    static public void printObjects(String label, Object[] objects)
    {
        if(objects == null) return;
        System.out.print(label);
        for(int n = 0; n < objects.length; n++)
            System.out.print(" " + objects[n].toString());
        System.out.println();
    }

    static public String toString(float[] values, String message)
    {
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append(message);
        if(values == null || values.length == 0)
            stringBuffer.append("no values");
        for(int n = 0; n < values.length; n++)
            stringBuffer.append(" " + values[n]);
        stringBuffer.append("\n");
        return stringBuffer.toString();
    }

    static public String toString(double[] values, String message)
    {
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append(message);
        if(values == null || values.length == 0)
            stringBuffer.append("no values");
        for(int n = 0; n < values.length; n++)
            stringBuffer.append(" " + values[n]);
        stringBuffer.append("\n");
        return stringBuffer.toString();
    }

    static public String toString(int[] values, String message)
    {
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append(message);
        if(values == null || values.length == 0)
            stringBuffer.append("no values");
        for(int n = 0; n < values.length; n++)
            stringBuffer.append(" " + values[n]);
        stringBuffer.append("\n");
        return stringBuffer.toString();
    }

    static public String getString(Object object)
    {
        if(object == null) return "null";
        return object.toString();
    }

    static public void main(String[] args)
    {
        Object object = new StsWiggleDisplayProperties();
        String simpleName = getSimpleClassname(object);
        System.out.println(simpleName);
        /*
            StsWiggleDisplayProperties wigglesDefaults = new StsWiggleDisplayProperties();

            StsWiggleDisplayProperties wiggles1 = new StsWiggleDisplayProperties(null, wigglesDefaults, "name");
            wiggles1.setOffsetAxisType(StsWiggleDisplayProperties.OFFSET_AXIS_VALUE);
            boolean isChanged = wiggles1.isChanged();
            System.out.println("isChanged " + isChanged);
        */
        /*
        TestCompareField field = new TestCompareField(1, 2.0f);
        TestCompare object1 = new TestCompare(3, 4.0f, field);
        TestCompare object2 = new TestCompare(3, 4.0f, field);
        boolean compares = compareObjects(object1, object2);
        System.out.println("Compares: " + compares);
        */
        /*
        beep();
              StsSeismicVolume volume = new StsSeismicVolume();
              StsVirtualVolume virtualVolume = new StsVirtualVolume();
              StsToolkit.copyDifferentClasses(volume, virtualVolume);
              System.out.println("done");
        */
        /*
                String packageName = args[0];
                String[] filenames = new String[] { "com/Sts/WorkflowPlugIn/PlugIns/StsWebSeis3dWorkflow.class" };
                ArrayList classNamesArray = truncateClassnames(filenames, "com.Sts.WorkflowPlugIn.PlugIns");
                if(classNamesArray == null || classNamesArray.size() == 0)
           System.out.println("No classes found.");
                else
                {
           for(int n = 0; n < classNamesArray.size(); n++)
           {
               String className = (String)classNamesArray.get(n);
               System.out.println("class: " + className);
           }
                }
        */
    }
}

class TestCompare
{
    int i;
    float f;
    TestCompareField fieldObject;

    TestCompare(int i, float f, TestCompareField fieldObject)
    {
        this.i = i;
        this.f = f;
        this.fieldObject = fieldObject;
    }
}

class TestCompareField
{
    int i;
    float f;

    TestCompareField(int i, float f)
    {
        this.i = i;
        this.f = f;
    }
}