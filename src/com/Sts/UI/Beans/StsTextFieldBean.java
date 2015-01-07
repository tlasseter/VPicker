package com.Sts.UI.Beans;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

abstract public class StsTextFieldBean extends StsFieldBean
{
    protected JTextField textField = new JTextField(nDefaultColumns);
	int nColumns = nDefaultColumns;
	static final int nDefaultColumns = 10;

    static final boolean debug = false;

    protected StsTextFieldBean()
    {
		setAlignment();
		addListeners();
    }
/*
    public StsTextFieldBean(boolean editable)
    {
        this(null, editable);
    }
*/
    protected StsTextFieldBean(boolean editable, String label)
    {
        super(editable, label);
		addListeners();
		setAlignment();
		setEditable();
    }

	public StsTextFieldBean(Object beanObject, String fieldName, String label)
	{
		initialize(beanObject, fieldName, true, label);
		addListeners();
		setAlignment();
		setEditable();
	}

	abstract protected void layoutBean();

/*
    public StsTextFieldBean(Class c, String fieldName, String label)
    {
        super();
        classInitialize(c, fieldName, true, label);
    }

    public StsTextFieldBean(Class c, String fieldName, boolean editable, String label)
    {
        super();
        classInitialize(c, fieldName, editable, label);
    }
*/
    protected void initialize(Object object, String fieldName, boolean editable, String label)
	{
		super.initialize(object, fieldName, editable, label);
		setEditable();
    }

    public void classInitialize(Class c, String fieldName, boolean editable, String label)
    {
        super.classInitialize(c, fieldName, editable, label);
		setEditable();
    }
/*
    public void classInitialize(String label)
    {
        layoutBean(false);
        setLabel(label);
    }
*/

    private void addListeners()
    {
        try
        {
            addActionListener();
            addFocusListener();
//            addKeyListener(); // checks validity on each key stroke
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    protected void addActionListener()
    {
        textField.addActionListener(this);
    }

    public void setToolTipText(String tip)
    {
        textField.setToolTipText(tip);
    }

	public void setBeanFont(Font font)
	{
		if(textField != null) textField.setFont(font);
        super.setBeanFont(font);
    }

    // override in any subclass with String values to turn off up/down edit arrows
    protected boolean isStringBean() { return false; }

    public void setPreferredSize(Dimension dimension)
    {
        textField.setPreferredSize(dimension);
    }

    public void setMinimumSize(Dimension dimension)
    {
        textField.setMinimumSize(dimension);
    }

    public void setMaximumSize(Dimension dimension)
    {
		textField.setMaximumSize(dimension);
	}

	public void setColumns(int nColumns)
	{
		this.nColumns = nColumns;
		textField.setColumns(nColumns);
//		textField.setHorizontalAlignment(JTextField.RIGHT);
	}
/*
    public void addListener(ActionListener listener)
    {
        if(actionListener != null) textField.removeActionListener(actionListener);
        actionListener = listener;
        textField.addActionListener(this);
        textField.setActionCommand(fieldName);
    }
*/
    public String getText() { return textField.getText(); }

    public void setText(String text)
    {
        final String t = text;
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { textField.setText(t); } } );
        // Runnable runnable = new Runnable() {  public void run()  {  textField.setText(t);  } };
        // StsToolkit.paintImmediately(runnable, textField);
    }

//    public Object getValueObject() { return textField.getText(); }

    public void actionPerformed(ActionEvent e)
    {
        boolean changed = false;
        if(e.getSource() != textField) return;

		boolean hasFocus = textField.hasFocus();
		String text = e.getActionCommand();
		if(!isValidString(text))
		{
			StsToolkit.beep();
//                textField.requestFocus();
		}
		else
		{
			setValueObject(text);
		}
        setValueInPanelObjects();
//        textField.setEnabled(false);
    }
/*
    protected void addKeyListener()
    {
        textField.addKeyListener
        (
            new KeyAdapter()
            {
                public void keyTyped(KeyEvent e)
                {
                    boolean hasFocus = textField.hasFocus();
                    char c = e.getKeyChar();
                    String newString = new String(getText() + c);
                    if(!isValidString(newString))
                    {
                        StsToolkit.beep();
                        e.consume();
                    }
                }
            }
        );
    }
*/
    protected void addFocusListener()
    {
        textField.addFocusListener
        (
            new FocusAdapter()
            {
                public void focusLost(FocusEvent e)
                {
                    Object source = e.getSource();
                    JTextField textField = (JTextField)e.getSource();
                    if(debug) System.out.println("focus lost for fieldBean " + label.getText());
                    boolean hasFocus = textField.hasFocus();
                    if(!isValidString(textField.getText()) )
                    {
                        StsToolkit.beep();
//                        textField.requestFocus();
                    }
                    else
                    {
                        String text = textField.getText();
                        setValueObject(text);
                        setValueInPanelObjects();
                        if(debug) System.out.println("fieldBean " + label.getText() + " value set to " + text);
                    }
                }
                public void focusGained(FocusEvent e)
                {
//                    textField.setEnabled(true);
                }
            }
        );
    }

    // Override these in subClasses
    public boolean decrement() { return false; }
    public boolean  increment() { return false; }

    public void doSetValueObject(Object textObject)
    {
        try
        {
            if(!(textObject instanceof String)) return;
            String text = (String)textObject;

            String oldText = textField.getText();
            boolean changed = !oldText.equals(text);
            if(!changed) return;
			setText(text);
			textField.setSelectionStart(0);
        }
        catch(Exception e)
        {
            StsToolkit.beep();
            textField.requestFocus();
        }
    }

	public void setAlignment()
	{
		textField.setScrollOffset(0);
	}

	public void setEditable()
	{
		textField.setEnabled(true);
		textField.setEditable(editable);
        if(editable)
        {
            textField.setBackground(Color.WHITE);
            textField.setBorder(BorderFactory.createLoweredBevelBorder());
//			addListeners();
        }
        else
        {
            textField.setBackground(Color.LIGHT_GRAY);
            textField.setBorder(BorderFactory.createEtchedBorder());
        }
    }

    public Object getValueObject() { return textField.getText(); }

    public Component[] getBeanComponents()
	{
		if(label != null)
			return new Component[] { label, textField };
		else
			return new Component[] { textField };
	}

	public Component getMainComponent() { return textField; }
    public JTextField getTextField() { return textField; }
    
    public boolean isValidChar(char c) { return true; }
    public boolean isValidString(String string) { return true; }
	public void print() { System.out.println("bean value " + textField.getText()); }
}

