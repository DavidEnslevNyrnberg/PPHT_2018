package com.example.erik.acutestressapp;



import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.view.View;
import android.widget.ArrayAdapter;
import android.widget.ListView;
import android.widget.TextView;

import java.util.ArrayList;
import java.util.List;


public class alarmHistory extends AppCompatActivity {

	// Creates the layout
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.alarm_history);
		final ListView listview = findViewById(R.id.listview);

		String[] values = new String[] { "08.04.18 - 12:25", "08.04.18 - 14:35", "08.04.18 - 16:10",
				"10.04.18 - 11:05", "01.05.18 - 22:10",   };
		ArrayAdapter<String> adapter = new ArrayAdapter<String>(this,R.layout.list_white_text,
				R.id.list_content, values);
		listview.setAdapter(adapter);


	}
}