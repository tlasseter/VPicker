package com.Sts.WorkflowPlugIn;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.text.*;
import javax.swing.text.html.*;
import java.awt.*;
import java.awt.event.*;
import java.net.*;
import java.util.*;

public class StsWorkflowNew extends JPanel
{
    protected StsWin3d win3d;
	protected JEditorPane editPane;
	protected HashMap m_controlMap = new HashMap();

	private StsWorkflowNew(URL editorPageURL, StsWin3d win3d) throws Exception
    {
        this.win3d = win3d;

		editPane = new JEditorPane();
		editPane.setEditable(false);
		editPane.setEditorKit(new HTMLEditorKit()
        {
			public ViewFactory getViewFactory()
            {
				return new TrackingFactory(m_controlMap);
			}
		});
        editPane.setPage(editorPageURL);
		JScrollPane scroll = new JScrollPane(editPane);
        setLayout(new BorderLayout());
		add(scroll, BorderLayout.CENTER);
	}

    static public StsWorkflowNew constructor(URL editorPageURL, StsWin3d win3d)
    {
        try
        {
            return new StsWorkflowNew(editorPageURL, win3d);
        }
        catch(Exception e)
        {
            if(editorPageURL != null)
                StsException.systemError("StsWorkflowEditor.constructor() failed for URL: " + editorPageURL.getPath());
            else
                StsException.systemError("StsWorkflowEditor.constructor() failed for URL: null");
            return null;
        }
    }

    private void submitWorkflow()
    {
        win3d.submitWorkflow();
    }

	private class TrackingFactory extends HTMLEditorKit.HTMLFactory
	{
		private HashMap m_map;

		public TrackingFactory(HashMap map)
        {
			m_map = map;
		}

		public View create(Element elem)
        {
			AttributeSet atts = elem.getAttributes();
			Object o = atts.getAttribute(StyleConstants.NameAttribute);
			String id = (String)atts.getAttribute(HTML.Attribute.ID);
			if(id != null && (HTML.Tag.INPUT == o || HTML.Tag.SELECT == o || HTML.Tag.TEXTAREA == o))
            {
				LinkedView view = new LinkedView(id, elem);
				m_map.put(id, view);
				return view;
			}
            else
            {
				return super.create(elem);
			}
		}
	}

	private class LinkedView extends FormView
	{
		private String m_id;

		public LinkedView(String id, Element elem)
        {
			super(elem);
			m_id = id;
		}

		protected Component createComponent()
        {
			Component comp = super.createComponent();
			if (comp instanceof JButton)
            {
				((AbstractButton)comp).addActionListener(new ActionListener()
                {
					public void actionPerformed(ActionEvent e)
                    {
						System.out.println("Event for button " + m_id);
                        if(m_id.equals("submit")) submitWorkflow();
					}
				});
			}
            else if (comp instanceof JComboBox)
            {
				((JComboBox)comp).addActionListener(new ActionListener()
                {
					public void actionPerformed(ActionEvent e)
                    {
						System.out.println("Selection field " + m_id + " is now " +
							((JComboBox)e.getSource()).getSelectedItem());
					}
				});
			}
            else if (comp instanceof JTextField)
            {
				((JTextField)comp).setInputVerifier(new InputVerifier()
                {
					public boolean verify(JComponent input)
                    {
						if (input instanceof JTextField)
                        {
							System.out.println("Text field " + m_id +
								" is now " + ((JTextField)input).getText());
						}
						return true;
					}
				});
			}
            else if (comp instanceof JCheckBox)
            {
				((JCheckBox)comp).addActionListener(new ActionListener()
                {
					public void actionPerformed(ActionEvent e)
                    {
						System.out.println("Checkbox " + m_id + " is checked " +
							((JCheckBox)e.getSource()).isSelected());
					}
				});
			}
			return comp;
		}
	}
}