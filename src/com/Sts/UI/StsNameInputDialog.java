package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsNameInputDialog extends JDialog
{
   public static boolean ALLOW_CANCEL = true;

   private StsJPanel panel = new StsJPanel();
   private JPanel buttonPanel = new JPanel();
   
   StsStringFieldBean nameBean = null;
   private JButton cancelButton = null;
   private JButton okayButton = new JButton("OK");

   private boolean wasCanceled = false;

   public StsNameInputDialog(Frame owner, String title, String inputLabel, String currentName, boolean allowCancel)
   {
      super(owner, title, true);
      try
      {
         setDefaultCloseOperation(DISPOSE_ON_CLOSE);

         nameBean = new StsStringFieldBean(this, "userName", inputLabel);
         nameBean.setValue(currentName);
         
         jbInit(allowCancel);
		 StsToolkit.centerComponentOnScreen(this);

         pack();
      }
      catch (Exception exception)
      {
         exception.printStackTrace();
      }
   }

   public String getUserName()
   {
      return nameBean.getValue();
   }

   public void setUserName(String name)
   {
      return;
   }
   public boolean wasCanceled()
   {
      return wasCanceled;
   }

   private void jbInit(boolean allowCancel)
      throws Exception
   {
      panel.addEndRow(nameBean,2,1.0);
      buttonPanel.setLayout(new FlowLayout());
      panel.addEndRow(buttonPanel);
      buttonPanel.add(okayButton);
      if (allowCancel)
      {
         cancelButton = new JButton("Cancel");
         buttonPanel.add(cancelButton);
         cancelButton.addActionListener(new java.awt.event.ActionListener()
         {
            public void actionPerformed(ActionEvent e)
            {
               cancelButtonAction();
            }
         });
      }
      okayButton.addActionListener(new java.awt.event.ActionListener()
      {
         public void actionPerformed(ActionEvent e)
         {
            okayButtonAction();
         }
      });

      getContentPane().add(panel);
   }

   private void okayButtonAction()
   {
      setVisible(false);
   }

   private void cancelButtonAction()
   {
      wasCanceled = true;
      setVisible(false);
   }

}
