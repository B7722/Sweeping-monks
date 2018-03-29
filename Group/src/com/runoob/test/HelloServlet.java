package com.runoob.test;

import java.io.IOException;
import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Servlet implementation class HelloServlet
 */
@WebServlet("/HelloServlet")
public class HelloServlet extends HttpServlet {
    private static final long serialVersionUID = 1L;
       
    /**
     * @see HttpServlet#HttpServlet()
     */
    public HelloServlet() {
        super();
        // TODO Auto-generated constructor stub
    }

    /**
     * @see HttpServlet#doGet(HttpServletRequest request, HttpServletResponse response)
     */
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        doPost(request,response);
    }

    /**
     * @see HttpServlet#doPost(HttpServletRequest request, HttpServletResponse response)
     */
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        // TODO Auto-generated method stub
    	String button = request.getParameter("button");
    	if(button.equals("Add")) {
	        String tetrahedron_1 = request.getParameter("tetrahedron_1");
	        String ball_p = request.getParameter("ball_p");
	        String cuboid_p = request.getParameter("cuboid_p");
	
	        String color_r = request.getParameter("color_r");
	        String color_g = request.getParameter("color_g");
	        String color_b = request.getParameter("color_b");
	        
	        String camera_p = request.getParameter("camera_p");
	        String camera_d = request.getParameter("camera_d");
	        String range = request.getParameter("range");
	
	        String txt = "Name: Camera"+"\r\n"+"Position: "+camera_p+"\r\n"+"Orientation: "+camera_d+"\r\n";
	        if(!tetrahedron_1.equals("")) {
	            String tetrahedron_2 = request.getParameter("tetrahedron_2");
	            String tetrahedron_3 = request.getParameter("tetrahedron_3");
	            String tetrahedron_4 = request.getParameter("tetrahedron_4");
	            txt = txt+"Name: Tetrahedran"+"\r\n"+"Point1: "+tetrahedron_1+"\r\n"+"Point2: "+tetrahedron_2+"\r\n"+"Point3: "+tetrahedron_3+"\r\n"+"Point4: "+tetrahedron_4+"\r\n"+"Color: "+color_r+" "+color_g+" "+color_b+"\r\n"+"Transparency: "+range+"\r\n";
	        }else if(!ball_p.equals("")) {
	        	String ball_r = request.getParameter("ball_r");
	        	txt = txt+"Name: Ball"+"\r\n"+"Position: "+ball_p+"\r\n"+"Radius: "+ball_r+"\r\n"+"Color: "+color_r+" "+color_g+" "+color_b+"\r\n"+"Transparency: "+range+"\r\n";
	        }else if(!cuboid_p.equals("")) {
	            String cuboid_x = request.getParameter("cuboid_x");
	            String cuboid_y = request.getParameter("cuboid_y");
	            String cuboid_z = request.getParameter("cuboid_z");
	        	txt = txt+"Name: Cuboid"+"\r\n"+"Position: "+cuboid_p+"\r\n"+"X: "+cuboid_x+"\r\n"+"Y: "+cuboid_y+"\r\n"+"Z: "+cuboid_z+"\r\n"+"Color: "+color_r+" "+color_g+" "+color_b+"\r\n"+"Transparency: "+range+"\r\n";
	        }
	        TxtFile.saveAsFileWriter(txt); 
	        request.getRequestDispatcher("menu.html").forward(request, response);//转发到成功页面
        } 
    	else if(button.equals("Clean")) {
    		DeleteTxt.delete();
	        request.getRequestDispatcher("menu.html").forward(request, response);//转发到成功页面
    	}
    	else if(button.equals("Save")){

	        request.getRequestDispatcher("menu.html").forward(request, response);//转发到成功页面
    	}
    }

}
