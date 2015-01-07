package com.Sts.UI.Beans;

import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.reflect.*;
import java.util.*;

/** A fieldBean provides a link between the member of a class and a GUI component which represents it.
 *  When constructed, the user can designate the class or a specific instance of the class.
 *  The user designates the class (rather than an instance) if in fact the instance is unknown and will
 *  be assigned later, or the instance will change over time.  The user designates the instance if the
 *  bean will always be tied to that instance (though in principal it could be changed).
 *
 *  A fieldBean is subclassed from StsJPanel and contains a component such as a JTextField or JComboBox and a label.
 *
 *  The user must provide a getter and optionally (but usually) a setter (if the bean is editable).
 *  The name of the getter and setter follows the usual convention with "get" or "set" prefixed to the name
 *  of the field (fieldName in the constructor) with the first letter of the fieldName capitalized.
 *
 *  A bean can also be associated with a beanPanel which in turn can be associated with one or more
 *  instances of the class the bean is associated with.  When the bean value is changed, it in turn changes the
 *  value of the instance or group of instances associated with the beanPanel.  This technique is used on
 *  The StsObjectPanel which is associated with the instances of a particular class which have been selected.
 *  Multiselecting a number of instances and changing the bean value on the StsObjectPanel changes the value in
 *  all the selected instances.
 *
 *  A beanPanel must implement the StsBeanPanelI interface which returns the panelObject or array of panelObjects
 *  which have been selected on the panel for the bean to operate on.
 *
 *  @see StsBeanPanelI
 * 
 */
abstract public class StsFieldBean extends StsJPanel implements ActionListener, TableCellRenderer, TableCellEditor, Cloneable
{
    /** panel this bean is optionally associated with */
    protected StsBeanPanelI beanPanel;
    /** instance this bean is optionally associated with */
    protected Object beanObject = null;
    /** class this bean is associated with. The beanPanel represents an instance or instances of this class; and beanObject is an instance of this class. */
    protected Class beanObjectClass = null;
    /** Name of the field in the beanObjectClass which this bean is associated with. */
    protected String fieldName = null;
    /** getMethod in the beanObjectClass for this bean. */
    protected Method get = null;
    /** setMethod in the beanObjectClass for this bean. */
    protected Method set = null;
    /** Not used anymore I don't think.  Getter is always required. */
    protected Field field = null;
    // True if field is editable using the bean in which case getter and setter must be defined.
    // False if field is not editable, but value is isVisible; getter must be defined. */
    protected boolean editable = false;
    /** optional label for the component on this bean */
    protected JLabel label = null;
    /** Should remove.  ActionListener for most components is this and is assigned by the fieldBean to the component.
     *  Other beans use ListSelectionListeners, etc.
     */
    protected ActionListener actionListener;
    /** GridBagConstraints grid width for this fieldBean on an StsJPanel; override if the bean should occupy more columns. */
    protected int gridwidth = 1;

	public StsFieldBean()
	{
		super(false);
	}

    public StsFieldBean(boolean editable, String label)
	{
		super(false);
		this.editable = editable;
		setLabel(label);
	}

	public StsFieldBean copy(Object beanObject)
	{
		StsException.systemError(this, "copy", "Developer: you need to implement this method in this object!");
		return null;
	}

	protected void layoutBean() { }

	protected void classInitialize(Class c, String fieldName, String fieldLabel)
	{
		this.classInitialize(c, fieldName, true, fieldLabel);
	}

    protected void classInitialize(Class c, String fieldName, boolean editable, String fieldLabel)
	{
		this.editable = editable;
        beanObjectClass = c;
        setLabel(fieldLabel);
        enableComponents();
		setEditable();
		this.fieldName = fieldName;
		if (fieldName == null || c == null)return;

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
                    StsException.systemError(this, "classInitialize", "The field " + fieldName + " in class " + StsToolkit.getSimpleClassname(c) + " has no get accessor or doesn't exist.");
					editable = false;
				}
			}
			Class returnType = null;
			if(get != null)
				returnType = get.getReturnType();
			else if(field != null)
				returnType = field.getType();

            set = StsToolkit.getAccessor(c, fieldName, "set", returnType);
			if (set == null && editable)
			{
				System.err.println("The field " + fieldName + " in class " + StsToolkit.getSimpleClassname(c) + " has no set accessor or doesn't exist.");
				editable = false;
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	protected void initialize(Object beanObject, String fieldName)
	{
		this.initialize(beanObject, fieldName, true, null);
	}

	protected void initialize(Object beanObject, String fieldName, boolean editable, String fieldLabel)
	{
		if (beanObject == null)
		{
			StsException.systemError("StsFieldBean.classInitialize() failed. beanObject cannot be null. Field name = " + fieldName);
			return;
		}
		this.editable = editable;
		enableComponents();
		setEditable();
		this.beanObject = beanObject;
		beanObjectClass =  beanObject.getClass();
		classInitialize(beanObjectClass, fieldName, editable, fieldLabel);
		setValueFromPanelObject(beanObject);
	}

	public void setEditable(boolean editable)
	{
		if (this.editable == editable)return;
		this.editable = editable;               
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { setEditable(); } } );
	}

    public boolean getEditable() { return editable; }

	public void setLabel(String fieldLabel)
	{
		if (fieldLabel == null) 
			return;
		label = new JLabel(fieldLabel);
		label.setHorizontalAlignment(SwingConstants.LEFT);
	}

	public String getLabelText()
	{
		if(label == null) return null;
		return label.getText();
	}

	public void resetLabel(String fieldLabel)
	{
		setLabel(fieldLabel);
		getParent().repaint();
	}

	protected void labelAndComponentLayout(Component label, Component component)
	{
		if (label != null)
		{
			gbc.anchor = GridBagConstraints.WEST;
			addToRow(label);
//			gbc.gridwidth = GridBagConstraints.REMAINDER;
		}
		gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.EAST;
		addToRow(component);
		gbc.anchor = GridBagConstraints.CENTER;
//		gbc.weightx = 1.0;
		gbc.fill = GridBagConstraints.NONE;
	}

	/*
	 public void labelAndComponentLayout(Component label, Component component)
	 {
	  if(label != null)
	  add(label, BorderLayout.WEST);
	  add(component, BorderLayout.EAST);
	 }

	 public void labelAndComponentLayout(Component label, Component component)
	 {
	  gbc.anchor = GridBagConstraints.WEST;
	  if(label != null)
	  {
	   gbc.weightx = 0.0;
	   gbc.fill = GridBagConstraints.NONE;
	   addToRow(label);
	   gbc.fill = GridBagConstraints.HORIZONTAL;
	   gbc.weightx = 1.0;
	  }
	  addToRow(component);
	  gbc.anchor = GridBagConstraints.CENTER;
	 }
	 */
	public StsBeanPanelI getBeanPanel()
	{
		return beanPanel;
	}

	public void setBeanPanel(StsBeanPanelI beanPanel)
	{
		this.beanPanel = beanPanel;
		getValueFromPanelObject();
	}

	public Object getBeanObject()
	{
		return beanObject;
	}

    /** For a bean whose class is defined but not it's instance, we can set the instance of for this bean from a beanObject which contains the member assigned to this bean.
     *
     * @param beanObject object containing the member corresponding to this bean (beanObject.member is instanace 
     */
    public void setBeanObject(Object beanObject)
	{
        if(beanObjectClass == null) return;
        // if we don't have a beanObject for this fieldBean and this call doesn't set it: complain!
        if(!beanObjectClass.isInstance(beanObject))
        {
            if(this.beanObject == null) // dont currently have a beanObject and probably wish that this one was being assigned here.
                System.out.println("this object of class " + StsToolkit.getSimpleClassname(beanObject) + " is not an instance of fieldBean assigned class " + beanObjectClass.getName());
            return;
        }
        this.beanObject = beanObject;
		setValueFromPanelObject(beanObject);
	}

    static public void setBeanObject(StsFieldBean[] fieldBeans, Object beanObject)
    {
        if(fieldBeans == null || beanObject == null) return;
        for(int n = 0; n < fieldBeans.length; n++)
            fieldBeans[n].setBeanObject(beanObject);
    }

    public void initializeBean()
    {
        if(beanObject == null) return;
        setValueFromPanelObject(beanObject);
    }

    public void setDefaultValueInPanelObject()
	{
		Object defaultObject = getDefaultValueObject();
		setValueObject(defaultObject);
		setValueInPanelObject();
    }

    public void setBeanObjectValue()
    {
        setValueInPanelObject();
    }
	public void setValueInPanelObject()
	{
		Object panelObject = getPanelObject();
		if (panelObject == null)return;
		Object valueObject = null;

		try
		{
			valueObject = getValueObject();
			if(valueObject == null) return;
			if (set != null)
			{
				Object[] valueObjects = new Object[] { valueObject };
				set.invoke(panelObject, valueObjects);
			}
			else if (field != null)
			{
				field.set(panelObject, valueObject);
				field = null;
			}
		}
		catch (Exception e)
		{
			outputFieldException(e, panelObject, valueObject);
		}
	}

    public void setValueInPanelObject(Object panelObject, Object valueObject)
    {
        if(field == null) return;
        try
        {
            field.set(panelObject, valueObject);
        }
		catch (Exception e)
		{
			outputFieldException(e, panelObject, valueObject);
		}
    }

    protected void outputFieldException(Exception e, Object panelObject, Object valueObject)
	{
		String objectString, valueString;
		if (panelObject == null)
			objectString = " for object: null";
		else
			objectString = " for object: " + StsToolkit.getSimpleClassname(panelObject);

		if (valueObject == null)
			valueString = " value: null";
		else
			valueString = " value: " + StsToolkit.getSimpleClassname(valueObject);

		if (set != null)
		{
			StsException.outputWarningException(this, "setValueInPanelObject", "setter " + set.getName() + objectString + valueString, e);
//            set = null;
		}
		else if (field != null)
		{
			StsException.systemError("StsFieldBean.setValueInPanelObject() failed.\n" + " setter not defined for field " + field.getName() + objectString + valueString);
//            field = null;
		}
		else
			StsException.systemError("StsFieldBean.setValueInPanelObject() failed.\n" + " no setter or field defined " + objectString + valueString);
	}

	public void setValueInPanelObjects()
	{
		Object[] panelObjects = getPanelObjects();
		if (panelObjects == null)return;
		Object panelObject = null;
		Object value = null;

		try
		{
			if (set != null)
			{
				value = getValueObject();
				if (value == null)return;
				//Object[] values = new Object[] {value};
				for (int n = 0; n < panelObjects.length; n++)
				{
					panelObject = panelObjects[n];
					if (panelObject != null) set.invoke(panelObject, value);
				}
			}
			else if (field != null)
			{
				value = getValueObject();
				if (value == null)return;
				for (int n = 0; n < panelObjects.length; n++)
				{
					panelObject = panelObjects[n];
					if (panelObject != null) field.set(panelObject, value);
				}
			}
		}
		catch (Exception e)
		{
			outputFieldException(e, panelObject, value);
		}
	}

	public void getValueFromPanelObject()
	{
		Object panelObject = getPanelObject();
		if (panelObject == null)return;
		setValueFromPanelObject(panelObject);
	}

	public void setValueFromPanelObject(Object panelObject)
	{
		Object valueObject = null;

		if (field == null && get == null)return;

		try
		{
            beanObject = panelObject;
            if (field != null) valueObject = field.get(panelObject);
			else valueObject = get.invoke(panelObject, new Object[0]);

			if (valueObject != null)
				setValueObject(valueObject);
			else
				this.setValueInPanelObject();
		}
		catch (Exception e)
		{
			if (valueObject != null)
				StsException.systemError(this, "setValueFromPanelObject", " failed for valueObject: " + StsToolkit.getSimpleClassname(valueObject) + e.getMessage());
			else
				StsException.systemError(this, "setValueFromPanelObject", " failed for valueObject: null " + e.getMessage());
		}
	}

    public Object getValueFromPanelObject(Object panelObject)
     {
         if (field == null && get == null) return null;

         try
         {
             if (field != null)
                 return field.get(panelObject);
             else
                 return get.invoke(panelObject, new Object[0]);
         }
         catch (Exception e)
		{
		    StsException.systemError(this, "getValueFromPanelObject", " field " + field.getName());
            return null;
        }
     }

    public Object getPanelObject()
	{
		if (beanPanel != null)
			return beanPanel.getPanelObject();
		else
			return beanObject;
	}

	public Object[] getPanelObjects()
	{
		if (beanPanel != null)
			return beanPanel.getPanelObjects();
        else if (beanObject != null)
			return new Object[] {beanObject};
		else
			return null;
	}

	public void enableComponents()
	{
		Component[] components = getBeanComponents();
		components = removeNullComponents(components);
		for (int n = 0; n < components.length; n++)
			if(components[n] != label) components[n].setEnabled(editable);
	}

	public String getFieldName()
	{
		return fieldName;
	}

	public String getLabelString()
	{
        if(label == null) return "";
		return label.getText();
	}

	public void updateLabel(String text)
	{
		label.setText(text);
	}

	public JLabel getLabel()
	{
		return label;
	}

	abstract public Object getValueObject();
	public Object getDefaultValueObject() { return null; } // override in subclasses as needed

	/** converts String to valueObject of this class.  Used in retreiving from user preferences file. */
	abstract public Object fromString(String string);

	/** converts current valueObject of this class to a string.  Used in saving to user preferences file. */
	abstract public String toString();

	public void setValueObjectFromString(String string)
	{
		Object valueObject = fromString(string);
		this.setValueObject(valueObject);
        this.setValueInPanelObject();
	}

	// override this in subclasses if you want changes in behaviour/looks between editable/noneditable
	protected void setEditable() { }

    public boolean setValueObject(Object valueObject_)
    {
        final Object valueObject = valueObject_;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                doSetValueObject(valueObject);
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
        return true;
    }

    public abstract void doSetValueObject(Object valueObject);

    public abstract Component[] getBeanComponents(); // gets label and main component
	public Component getMainComponent() { return this; } // override in subclasses as needed

	public Component[] getNonNullComponents()
	{
		Component[] components = getBeanComponents();
		components = removeNullComponents(components);
		return components;
	}

	private Component[] removeNullComponents(Component[] components)
	{
		if (components == null)return new Component[0];
		int nNotNull = 0;
		for (int n = 0; n < components.length; n++)
			if (components[n] != null) nNotNull++;

		if (nNotNull == components.length)return components;

		Component[] newComponents = new Component[nNotNull];
		int nn = 0;
		for (int n = 0; n < components.length; n++)
			if (components[n] != null) newComponents[nn++] = components[n];
		return newComponents;
	}

	protected void addActionListener()
	{
		actionListener = this;
	}

    public void setBeanFont(Font font)
	{
		if (label != null) label.setFont(font);
    }

	// override these in subclass as needed
	public void actionPerformed(ActionEvent e)
	{
		StsException.systemError(this, "actionPerformed", "Need to implement action performed in subclass");
	}

	static public StsFieldBean getBeanWithFieldName(StsFieldBean[] beans, String name)
	{
		if (beans == null)
		{
			StsException.systemError("Unable to locate bean with name: " + name + " beans are null.");
			return null;
		}
		for (int n = 0; n < beans.length; n++)
			if (beans[n].getFieldName().equals(name))return beans[n];

		StsException.systemError("Unable to locate bean with name: " + name + " beans are null.");
		return null;
	}

    static public boolean getBeanWithFieldNameValue(StsFieldBean[] beans, String name, Object value)
    {
        StsFieldBean fieldBean = getBeanWithFieldName(beans, name);
        if(fieldBean == null) return false;
        fieldBean.setValueObject(value);
        return true;
    }

    public void outputSetValueObjectException(Object object, String correctValueObjectClassname)
	{
		StsException.systemError(this, "setValueObject",  "Object should be instance of " + correctValueObjectClassname +
								 " but instead it is an instance of " + object.getClass().getName());
	}

	public void addToPanel(StsFieldBeanPanel beanPanel)
	{
		GridBagConstraints gbc = beanPanel.gbc;
		gbc.gridwidth = gridwidth;
		Component[] beanComponents = getBeanComponents();
		gbc.gridx = 0;
		for (int n = 0; n < beanComponents.length; n++)
			beanPanel.addToRow(beanComponents[n]);
		gbc.gridy++;
		gbc.gridwidth = 1;
	}

	public void addComponentsToPanel(StsJPanel panel)
	{
		Component component = getMainComponent();
		GridBagConstraints gbc = panel.gbc;
		int gbcFill = gbc.fill; // save fill for restoring at end
		int gbcAnchor = gbc.anchor;
		gbc.gridwidth = gridwidth;
		gbc.fill = GridBagConstraints.NONE;
		gbc.weightx = 0.0;
		gbc.anchor = GridBagConstraints.WEST;
		if(label != null)
			panel.addToRow(label);
		gbc.fill = gbcFill;
		gbc.weightx = 1.0;
		if(gbcFill == GridBagConstraints.NONE)
			gbc.fill = GridBagConstraints.BOTH;
		if(component != null)
			panel.addToRow(component);
		gbc.fill = gbcFill;
		gbc.anchor = gbcAnchor;
	 }

	public void addToPanel(StsObjectPanel objectPanel)
	{
		objectPanel.add(this);
	}

	public Component getObjectPanelComponent()
	{
		StsJPanel panel = StsJPanel.addInsets();
		panel.objectPanelLayout(getBeanComponents());
		return panel;
	}

	public String getBeanKey()
	{
		if(beanObject == null) return null;
		return beanObject.getClass().getName() + "." + fieldName;
	}

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,  int row, int column)
    {
        setValueObject(value);
        return getMainComponent();
    }

    public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column)
    {
        setValueObject(value);
        return getMainComponent();
    }

    public Object getCellEditorValue()
    {
        return getValueObject();
    }

    public boolean isCellEditable(EventObject anEvent)
    {
        return editable;
    }

    public boolean shouldSelectCell(EventObject anEvent)
    {
        return editable;
    }

    public boolean stopCellEditing()
    {
        return true;
    }

    public void cancelCellEditing()
    {
    }

    public void addCellEditorListener(CellEditorListener l)
    {
    }

    public void removeCellEditorListener(CellEditorListener l)
    {        
    }
}
