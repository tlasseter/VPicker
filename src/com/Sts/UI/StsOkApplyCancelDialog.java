package com.Sts.UI;

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/** This is a dialog centered on a parent with a displayPanel displayed and an okCancelObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */

public class StsOkApplyCancelDialog extends JDialog
{
	public StsJPanel panel = StsJPanel.addInsets();
	public StsDialogFace[] okCancelObjects;
	public StsGroupBox buttonBox = new StsGroupBox();
	public StsButton okButton = new StsButton("OK", "Accept changes and dismiss dialog.", this, "ok");
	public StsButton applyButton = new StsButton("Apply", "Apply changes and keep dialog.", this, "apply");
	public StsButton cancelButton = new StsButton("Cancel", "Reject changes and dismiss dialog.", this, "cancel");
	private Object[] buttonRowObjects = new Object[] { okButton, applyButton, cancelButton };

    public StsOkApplyCancelDialog(Frame frame, String title, boolean modal)
    {
        super(frame, title, modal);
    }

	public StsOkApplyCancelDialog(Frame frame, StsDialogFace[] okCancelObjects, String title, boolean modal)
	{
		this(frame, title, modal);
		if(frame != null) setLocation(frame.getLocation());
		this.okCancelObjects = okCancelObjects;
		layoutPanels(false);
		addWindowCloseOperation();
		display();
	}

	public StsOkApplyCancelDialog(Frame frame, StsDialogFace[] okCancelObjects, String title, boolean modal, boolean expandIt)
	{
		this(frame, title, modal);
		if(frame != null) setLocation(frame.getLocation());
		this.okCancelObjects = okCancelObjects;
		layoutPanels(expandIt);
		addWindowCloseOperation();
		display();
	}

    public StsOkApplyCancelDialog(Frame frame, StsDialogFace okCancelObject, String title, boolean modal)
	{
		this(frame, new StsDialogFace[] { okCancelObject }, title, modal);
	}

    public void display()
	{    	
		pack();					
		setVisible(true);	
	}

	private void layoutPanels(boolean expandIt)
	{
		panel.gbc.fill = GridBagConstraints.HORIZONTAL;
		for(int n = 0; n < okCancelObjects.length; n++)
        {
            okCancelObjects[n] = okCancelObjects[n].getEditableCopy();
            panel.add(okCancelObjects[n].getPanel(expandIt));
        }	
        Insets insets = new Insets(4, 4, 4, 4);

        okButton.setMargin(insets);
		cancelButton.setMargin(insets);
		applyButton.setMargin(insets);
		buttonBox.addToRow(okButton);
		buttonBox.addToRow(applyButton);
		buttonBox.addEndRow(cancelButton);

		panel.add(buttonBox);
		getContentPane().add(panel);
	}

    public void collapseAll()
    {
        for(int i=0; i<panel.getComponentCount(); i++)
        {
            if(panel.getComponent(i) instanceof StsPropertiesPanel)
                ((StsPropertiesPanel)panel.getComponent(i)).dismissPanel();
        }
    }

	public void addWindowCloseOperation()
	{
		setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		addWindowListener
		 (
			 new WindowAdapter()
			 {
				 public void windowClosing(WindowEvent e)
				 {
					 close();
			     }
			 }
		);
	}

	public void close()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.CLOSING);
		dispose();
	}

	public void ok()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.OK);
		dispose();
	}

    public void submit()
    {
        if(StsYesNoDialog.questionValue(this,"This will take some time. Do you wish to continue?"))
        {
            // Save all the current settings
            for(int n = 0; n < okCancelObjects.length; n++)
                okCancelObjects[n].dialogSelectionType(StsDialogFace.OK);
            // Dispose of the dialog
            dispose();
            // Begin the processing

        }
	}

	public void cancel()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.CANCEL);
		setVisible(false);
		dispose();
	}

	public void apply()
	{
		for(int n = 0; n < okCancelObjects.length; n++)
			okCancelObjects[n].dialogSelectionType(StsDialogFace.APPLY);
    }


	static public void main(String[] args)
	{
		TestOkApplyCancel testObject = new TestOkApplyCancel();
		StsOkApplyCancelDialog dialog = new StsOkApplyCancelDialog((Frame)null, testObject, "Title.", true);
	}
}

class TestOkApplyCancel implements StsDialogFace
{
	String string = "test";
	StsStringFieldBean stringBean = new StsStringFieldBean(this, "string", "String");

	public TestOkApplyCancel()
	{
	}

    public StsDialogFace getEditableCopy()
    {
        return new TestOkApplyCancel();
    }

    public void dialogSelectionType(int type)
	{
		System.out.println("Selection Type " + type);
	}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		StsGroupBox groupBox = new StsGroupBox("Test OkApplyCancel");
		groupBox.add(stringBean);
		return (Component)groupBox;
	}

	public void setString(String s) { string = s; }
	public String getString() { return string; }
}
