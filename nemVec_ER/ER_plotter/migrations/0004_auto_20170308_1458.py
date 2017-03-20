# -*- coding: utf-8 -*-
# Generated by Django 1.10.5 on 2017-03-08 13:58
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ER_plotter', '0003_fasta'),
    ]

    operations = [
        migrations.CreateModel(
            name='Cpm',
            fields=[
                ('nvertx_id', models.CharField(max_length=20, primary_key=True, serialize=False)),
                ('regen_anc_UC', models.FloatField(max_length=10)),
                ('regen_anc_0HPA', models.FloatField(max_length=10)),
                ('regen_anc_2HPA', models.FloatField(max_length=10)),
                ('regen_anc_4HPA', models.FloatField(max_length=10)),
                ('regen_anc_8HPA', models.FloatField(max_length=10)),
                ('regen_anc_12HPA', models.FloatField(max_length=10)),
                ('regen_anc_16HPA', models.FloatField(max_length=10)),
                ('regen_anc_20HPA', models.FloatField(max_length=10)),
                ('regen_anc_24HPA', models.FloatField(max_length=10)),
                ('regen_anc_36HPA', models.FloatField(max_length=10)),
                ('regen_anc_48HPA', models.FloatField(max_length=10)),
                ('regen_anc_60HPA', models.FloatField(max_length=10)),
                ('regen_anc_72HPA', models.FloatField(max_length=10)),
                ('regen_anc_96HPA', models.FloatField(max_length=10)),
                ('regen_anc_120HPA', models.FloatField(max_length=10)),
                ('regen_anc_144HPA', models.FloatField(max_length=10)),
            ],
        ),
        migrations.RemoveField(
            model_name='fasta',
            name='id',
        ),
        migrations.AlterField(
            model_name='fasta',
            name='nvertx_id',
            field=models.CharField(max_length=20, primary_key=True, serialize=False),
        ),
    ]
