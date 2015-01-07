package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.UI.DataTransfer.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p/>
 * <p>Description: </p>
 * <p/>
 * <p>Copyright: Copyright (c) 2006</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

abstract public class StsSimpleOkCancelDialog extends JDialog
{
    public static boolean ALLOW_CANCEL = true;
    private StsJPanel backgroundPanel = new StsJPanel();
    private StsJPanel panel = StsJPanel.addInsets();
    private StsJPanel buttonPanel = new StsJPanel();
    private StsButton cancelButton = null;
    private StsButton okayButton = new StsButton("OK", "Accept panel edits.", this, "ok");

    private boolean wasCanceled = false;

    abstract public void layoutPanel(StsJPanel panel);

    public StsSimpleOkCancelDialog(Frame owner, String title)
    {
        super(owner, title, true);

        try
        {
            setDefaultCloseOperation(DISPOSE_ON_CLOSE);
            constructPanel();

//            StsToolkit.centerComponentOnScreen(this);
            pack();
            setLocation(owner);
        }
        catch(Exception exception)
        {
            exception.printStackTrace();
        }
    }

    public StsSimpleOkCancelDialog(Frame owner, String title, int x, int y)
    {
        this(owner, title);

    }

    public void setLocation(Frame owner)
    {
        Point location = owner.getLocationOnScreen();
        Point offset = owner.getMousePosition();
        location.translate(offset.x, offset.y);
        super.setLocation(location);
    }

    public boolean wasCanceled()
    {
        return wasCanceled;
    }

    private void constructPanel() throws Exception
    {
        layoutPanel(panel);
        buttonPanel.addToRow(okayButton);
        cancelButton = new StsButton("Cancel", "Cancel: no changes.", this, "cancel");
        buttonPanel.addToRow(cancelButton);
        backgroundPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        backgroundPanel.add(panel);
        backgroundPanel.add(buttonPanel);
        getContentPane().add(backgroundPanel);
    }

    public void ok()
    {
        setVisible(false);
    }

    public void cancel()
    {
        wasCanceled = true;
        setVisible(false);
    }

    public static void main(String[] args)
    {
        TestSimpleOKCancelDialog d = new TestSimpleOKCancelDialog(null, "Test Dialog");
        d.setVisible(true);
        System.exit(0);
    }
}

class TestObjectsTransferPanel extends StsObjectsTransferPanel implements StsObjectTransferListener
{
    String[] objects = new String[]{"view one", "view two", "view three"};

    TestObjectsTransferPanel(String title, int width, int height)
    {
        initialize(title, this, width, height);
    }

    public Object[] getAvailableObjects() { return objects; }
    public Object[] initializeAvailableObjects() { return objects; }

    public void addObjects(Object[] objects)
    {
        System.out.println("add Objects:");
        printObjects(objects);
    }

    private void printObjects(Object[] objects)
    {
        for (int n = 0; n < objects.length; n++)
            System.out.println("    " + objects[n].toString());
    }

    public void removeObjects(Object[] objects)
    {
        System.out.println("remove Objects:");
        printObjects(objects);
    }

    public void objectSelected(Object selectedObject)
    {
        System.out.println("selected Object:" + selectedObject.toString());
    }
}

class TestSimpleOKCancelDialog extends StsSimpleOkCancelDialog
{
    TestObjectsTransferPanel transferPanel;

    TestSimpleOKCancelDialog(Frame frame, String title)
    {
        super(frame, title);
        Object[] selectedObjects = transferPanel.getSelectedObjects();
        System.out.print("selected objects:");
        for(int n = 0; n < selectedObjects.length; n++)
            System.out.print(" " + selectedObjects[n].toString());
        System.out.println();
    }

    public void layoutPanel(StsJPanel panel)
    {
        transferPanel = new TestObjectsTransferPanel("Test Panel", 400, 100);
        panel.add(transferPanel);
    }
}
