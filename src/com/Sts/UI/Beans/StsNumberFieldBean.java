package com.Sts.UI.Beans;

import com.Sts.Utilities.*;

import java.awt.*;
import java.awt.event.*;
import java.beans.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
abstract public class StsNumberFieldBean extends StsTextFieldBean
{
	protected ArrowBean up = null, down = null;
	boolean useArrows = false;
	protected StsJPanel arrowTextField;
	protected boolean continuousMode = true;

	public StsNumberFieldBean()
	{
	}

	public StsNumberFieldBean(boolean editable)
	{
		this(editable, null);
	}

	public StsNumberFieldBean(boolean editable, String label)
	{
		super(editable, label);
	}

	public StsNumberFieldBean(Object beanObject, String fieldName, String label)
	{
		super(beanObject, fieldName, label);
	}

	protected void layoutBean(boolean useArrows)
	{
		this.useArrows = useArrows;
		layoutBean();
	}

	protected void layoutBean()
	 {
		 try
		 {
			 if(editable && useArrows)
			 {
				 addArrows();
				 labelAndComponentLayout(label, arrowTextField);
			 }
			 else
				labelAndComponentLayout(label, textField);

//			 textField.setColumns(10);
//			 textField.setSize(50, 20);
		 }
		 catch (Exception e)
		 {
			 e.printStackTrace();
		 }
	 }

	public void setContinuousMode(boolean mode)
	{
		continuousMode = mode;
	}
	 private void addArrows()
	 {
		 try
		 {
			 arrowTextField = StsJPanel.addInsets();
			 GridBagConstraints gbc = arrowTextField.gbc;
			 // don't stretch up/down arrows horizontally
			 gbc.fill = GridBagConstraints.NONE;
			 gbc.anchor = GridBagConstraints.EAST;
			 gbc.weightx = 0.0;
			 down = (ArrowBean) Beans.instantiate (getClass().getClassLoader (), "com.Sts.UI.Beans.ArrowBean");
			 arrowTextField.addToRow((Component) Beans.getInstanceOf (down, Component.class));
			 down.setDirection (ArrowBean.LEFT);
			 down.setSize(10, 20);
			 down.addActionListener (this);

			 gbc.fill = GridBagConstraints.BOTH;
			 gbc.anchor = GridBagConstraints.CENTER;
			 gbc.weightx = 1.0;
			 arrowTextField.addToRow(textField);

			 gbc.fill = GridBagConstraints.NONE;
			 gbc.anchor = GridBagConstraints.WEST;
			 gbc.weightx = 0.0;
			 up = (ArrowBean) Beans.instantiate (getClass ().getClassLoader (), "com.Sts.UI.Beans.ArrowBean");
			 arrowTextField.addEndRow((Component) Beans.getInstanceOf (up, Component.class));
			 up.setDirection (ArrowBean.RIGHT);
			 up.setSize(10, 20);
			 up.addActionListener (this);
			 gbc.weightx = 1.0;
		}
		 catch (Exception e)
		 {
			 e.printStackTrace();
			 return;
		 }
    }
	public Component[] getBeanComponents()
	{
		Component mainComponent = getMainComponent();
		if(this.label != null)
			return new Component[] { label, mainComponent };
		else
			return new Component[] { mainComponent };
	}

	public Component getMainComponent()
	{
		if(this.arrowTextField != null)
			return arrowTextField;
		else
			return textField;
    }

	public void actionPerformed(ActionEvent e)
	{
//		 setEnabled(true);
         if(!getEditable())
             return;
		 boolean changed = false;
		 if (e.getSource () == down)
		 {
			 if(!continuousMode && down.down)
				 changed = false;
			 else
				 changed = decrement();
		 }
		 else if (e.getSource () == up)
		 {
			 if(!continuousMode && up.down)
				 changed = false;
			 else			 
				 changed = increment();
		 }
		 else if(e.getSource() == textField)
		 {
			 boolean hasFocus = textField.hasFocus();
			 String text = e.getActionCommand();
			 if(!isValidString(text))
			 {
				 StsToolkit.beep();
//                textField.requestFocus();
			 }
			 else
			 {
				// jbw 10/16/09 this prevents check range from firing doSetValueObject(text);
				setValueObject(text);
			 }
		 }
	     setValueInPanelObjects();
    }
}
