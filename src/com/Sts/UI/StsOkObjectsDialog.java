package com.Sts.UI;

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/** This is a dialog centered on a parent with a displayPanel displayed and an okObject which it communicates with
 *  when ok/apply/cancel buttons are pushed.
 */


public class StsOkObjectsDialog extends JDialog
{
    StsJPanel backPanel = StsJPanel.addInsets();
    public StsDialogFace[] okObjects;
	private StsGroupBox buttonBox = new StsGroupBox();
	private StsButton okButton = new StsButton("OK", "Accept changes and dismiss dialog.", this, "ok");

	public StsOkObjectsDialog(Frame frame, StsDialogFace[] okObjects, String title, boolean modal)
	{
		super(frame, title, modal);
        if(frame != null) setLocation(frame.getLocation());
		this.okObjects = okObjects;
        constructPanel(okObjects);
        getContentPane().add(backPanel);
		addWindowCloseOperation();
		display();
	}

	public StsOkObjectsDialog(Frame frame, StsDialogFace okCancelObject, String title, boolean modal)
	{
		this(frame, new StsDialogFace[] { okCancelObject }, title, modal);
	}

    public StsOkObjectsDialog(Frame frame, JPanel panel, String title, boolean modal)
	{
		super(frame, title, modal);
        if(frame != null) setLocation(frame.getLocation());
        constructPanel(panel);
		addWindowCloseOperation();
		display();
	}

    public void display()
	{
		pack();
		setVisible(true);
	}

    private StsJPanel constructPanel(StsDialogFace[] okObjects)
	{
		backPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		for(int n = 0; n < okObjects.length; n++)
			backPanel.add(okObjects[n].getPanel());
		Insets insets = new Insets(4, 4, 4, 4);
		okButton.setMargin(insets);
		buttonBox.addToRow(okButton);
		backPanel.add(buttonBox);
        return backPanel;
    }

    private void constructPanel(JPanel panel)
	{
		backPanel.add(panel);
		Insets insets = new Insets(4, 4, 4, 4);
		okButton.setMargin(insets);
		buttonBox.addToRow(okButton);
		panel.add(buttonBox);
		getContentPane().add(panel);
	}

	private void addWindowCloseOperation()
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
		for(int n = 0; n < okObjects.length; n++)
			okObjects[n].dialogSelectionType(StsDialogFace.CLOSING);
		dispose();
	}

	public void ok()
	{
        if(okObjects != null)
        {
            for(int n = 0; n < okObjects.length; n++)
                okObjects[n].dialogSelectionType(StsDialogFace.OK);
        }
        dispose();
	}

	static public void main(String[] args)
	{
		TestOkObj testObject = new TestOkObj();
		StsOkObjectsDialog dialog = new StsOkObjectsDialog(null, new StsDialogFace[] { testObject }, "Title.", true);
	}
}

class TestOkObj implements StsDialogFace
{
	String string = "test";
	StsStringFieldBean stringBean = new StsStringFieldBean(this, "string", "String");
	TestOkObj()
	{
	}

	public void dialogSelectionType(int type)
	{
		System.out.println("Selection Type " + type);
	}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		StsGroupBox groupBox = new StsGroupBox("Test OkCancel");
		groupBox.add(stringBean);
		return (Component)groupBox;
	}

    public StsDialogFace getEditableCopy()
    {
        return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
    }

    public void setString(String s) { string = s; }
	public String getString() { return string; }
}

