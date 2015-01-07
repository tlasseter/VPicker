
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Beans;

import com.Sts.Actions.Import.*;
import com.Sts.IO.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

public class StsParameterFileFieldBean extends StsFieldBean
{
    StsButton newButton = new StsButton();
    StsComboBoxFieldBean fileListBean = new StsComboBoxFieldBean();
    StsStringFieldBean newNameBean = new StsStringFieldBean();
    Class subClass = null;
    Class superClass = null;
    Object[] items;
    String btnTooltip = "Output a new parameter file.";
    String cboxTooltip = "Select an existing parameter file.";
    String name = "None";
    String format = null, group = null;
    StsJPanel ioObject = null;

    public StsParameterFileFieldBean()
    {
    }

    public StsParameterFileFieldBean(StsJPanel object, String group, String fmt)
    {
        this(object, object.getClass(), object.getClass(), group, fmt);
    }
    public StsParameterFileFieldBean(StsJPanel object, Class subClass, Class superClass, String group, String fmt)
    {
        this.subClass = subClass;
        this.superClass = superClass;
        this.ioObject = object;
        this.format = fmt;
        this.group = group;
        layoutBean();
    }

	public StsParameterFileFieldBean copy(Object beanObject)
	{
		return new StsParameterFileFieldBean(ioObject, subClass, superClass, group, format);
	}

    public void setFilename(String filename)
    {
        if(filename.equalsIgnoreCase("None"))
            return;
        this.name = filename;
        try
        {
            String name = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator +
                          group + "." + format + "." + this.name;
            StsParameterFile.initialReadObjectFields(name, ioObject, subClass, superClass);
            ioObject.setPanelObject(ioObject);
        }
        catch(Exception e)
        {
            Frame frame = new Frame("Read Error");
            frame.setLocation(ioObject.getLocation());
            new StsMessage(frame, StsMessage.ERROR, "Failed to read selected template: " + filename);
        }
    }

    public String getFilename()
    {
        return name;
    }

    public void setFormat(String fmt)
    {
        format = fmt;
        initialize();
    }
    public void setGroup(String grp)
    {
        group = grp;
        initialize();
    }

    public void newFilename()
    {
        // Remap bean to allow file name input
        if(newButton.getText().equals("New"))
        {
            this.remove(fileListBean);
            newNameBean.setValue("<Enter Name>");
            this.name = "None";
            newButton.setText("Save");
            add(newNameBean, BorderLayout.CENTER);
        }
        // Save file and remap bena to file selection
        else
        {
            // Save pressed but no name entered
            if(name.equals("None"))
                return;

            // Output parameters in named file
            if(outputFile(name))
            {
                this.remove(newNameBean);
                newButton.setText("New");
                add(fileListBean, BorderLayout.CENTER);
                initialize();
                fileListBean.setSelectedItem(name);
            }
            else
                newNameBean.setValue("<Enter Name>");
        }
        revalidate();
    }

    public boolean outputFile(String name)
    {
        // Output template
        try
        {
            String filename = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator  +
                          group + "." + format + "." + name;
            File file = new File(filename);
            if(file.exists())
            {
                if(!StsYesNoDialog.questionValue(this,"A file named " + name + " already exists. Overwrite it?"))
                    return false;
            }
            StsParameterFile.writeObjectFields(filename, ioObject, subClass, superClass);
            return true;
        }
        catch(Exception e)
        {
            Frame frame = new Frame("Write Error");
            frame.setLocation(ioObject.getLocation());
            if(StsYesNoDialog.questionValue(this,"Failed to write file name:" + name + ". Try Again?"))
                return false;
            else
                return true;

        }
    }

    public void initialize()
    {
        if(format == null || group == null)
            return;

        fileListBean.removeAll();

        String[] filenames = getFilenames();
        int nFiles = filenames.length;
        if (nFiles > 0)
        {
            fileListBean.addItem("None");
            for (int i = 0; i < filenames.length; i++)
                fileListBean.addItem(filenames[i]);
            fileListBean.setEditable(true);
        }
        else
        {
            fileListBean.addItem("None");
            fileListBean.setEditable(false);
        }
    }

    public String[] getFilenames()
    {
        String[] filenames;
        StsFilenameFilter filter = new StsFilenameFilter(group,format);

        // Create list of templates from Project directory
        String directory = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
        StsFileSet availableFileSet = StsFileSet.constructor(directory, filter);
        StsAbstractFile[] availableFiles = availableFileSet.getFiles();
        int nFiles = availableFiles.length;
        if (nFiles > 0)
        {
            filenames = new String[nFiles];
            for (int n = 0; n < nFiles; n++)
                filenames[n] = filter.getFilenameEnding(availableFiles[n].getFilename());
            return filenames;
        }
        else
            return new String[0];
    }

    protected void layoutBean()
    {
        try
        {
            fileListBean = new StsComboBoxFieldBean(this, "filename", "Available:");
            fileListBean.setToolTipText(cboxTooltip);
            newButton = new StsButton("New",btnTooltip,this,"newFilename");
            newNameBean = new StsStringFieldBean(this, "newName", "<Enter Name>", true, "New Name:");

            setLayout(new BorderLayout());
            add(fileListBean, BorderLayout.CENTER);
            add(newButton, BorderLayout.EAST);
            if (label != null)
                add(label, BorderLayout.WEST);
            initialize();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    public Component[] getBeanComponents()
    {
        return fileListBean.getBeanComponents();
    }

    public Object getValueObject()
    {
        return fileListBean.getComboBox().getSelectedItem();
    }

    public String toString()
    {
        return fileListBean.getComboBox().getSelectedItem().toString();
    }

    public Object fromString(String string)
    {
        if(items == null) return null;
        for(int n = 0; n < items.length; n++)
            if(items[n].toString().equals(string))
                return items[n];
        return null;
    }

    public String getNewName() { return "<Enter Name>"; }
    public void setNewName(String name)
    {
        this.name = name;
       // newFilename();
    }
    public void setListItems(Object[] items)
    {
        if (items == null || items.length == 0)
        {
            setEditable(false);
            return;
        }
        setEditable(true);
        this.items = items;
        DefaultComboBoxModel comboModel = new DefaultComboBoxModel(items);
        fileListBean.getComboBox().setModel(comboModel);
        comboModel.setSelectedItem(items[0]);
    }

    public void setSelectedItem(final Object item)
    {
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { fileListBean.getComboBox().setSelectedItem(item); }});
    }

    public void doSetValueObject(Object object)
    {
        Object item = getItem(object);
        if (item == null)return;
        fileListBean.getComboBox().getModel().setSelectedItem(item);
    }

    // override in subclass if object is inside an item
    public Object getItem(Object object)
    {
        if (items == null)return null;
        else return object;
    }

    public void setComboBoxToolTipText(String tip)
    {
        fileListBean.getComboBox().setToolTipText(tip);
    }

    public void setNewButtonToolTipText(String tip)
    {
        newButton.setToolTipText(tip);
    }

    public int getSelectedIndex()
    {
        return fileListBean.getComboBox().getSelectedIndex();
    }

    public void setSelectedIndex(int index)
    {
        fileListBean.getComboBox().setSelectedIndex(index);
    }

    public void setToLastItem()
    {
        int nItems = fileListBean.getComboBox().getModel().getSize();
        fileListBean.getComboBox().setSelectedIndex(nItems - 1);
    }

    public JComboBox getComboBox()
    {
        return fileListBean.getComboBox();
    }

    public void setPreferredSize(int width, int height)
    {
        fileListBean.getComboBox().setPreferredSize(new Dimension(width, height));
    }

    public void addItem(Object object)
    {
        setEditable(true);
        fileListBean.getComboBox().addItem(object);
    }

    protected void setEditable()
    {
        fileListBean.getComboBox().setEnabled(editable);
    }

    public void removeAll()
    {
        fileListBean.getComboBox().removeAllItems();
    }

    public void repaint()
    {
        if(fileListBean == null) return;
        if(fileListBean.getComboBox() != null)
            fileListBean.getComboBox().repaint();
    }

    static final class StsFilenameFilter implements FilenameFilter
    {
        String group;
        String format;
        int length;

        public StsFilenameFilter(String group, String format)
        {
            this.group = group;
            this.format = format;
            length = group.length() + format.length() + 2; // string length of "group.format."
        }

        public boolean accept(File dir, String name)
        {
            StsKeywordIO.parseAsciiFilename(name);
            String keywordGroup = StsKeywordIO.group;
            String keywordFormat = StsKeywordIO.format;
            return StsKeywordIO.group.equals(group) && StsKeywordIO.format.equals(format);
        }

        public String getFilenameEnding(String filename)
        {
            return filename.substring(length);
        }
    }
}
