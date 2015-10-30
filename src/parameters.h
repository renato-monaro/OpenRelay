# Copyright (C) 2008-2015 Renato Machado Monaro
# This file is part of OpenRelay.
#
# OpenRelay is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenRelay is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
#ifndef PARAMETERS_H
#define PARAMETERS_H

//! Parâmetro que define o tamanho do buffer circular.
//#define BUFFER_SIZE 512
//! Parâmetro que define qual núcleo do processador será utilizado.
#define CPUMAP 0x1
//! Parâmetro que define a aquisição com maior prioridade.
#define ACQUISITION_PRIO 0
//! Parâmetro que define a taxa de aquisição.
//#define ACQUISITION_TICK 520833  /*ns (1920Hz)*/
//! Parâmetro que define qual escala (range) será utilizada.
#define COMEDI_RANGE 0
//! Parâmetro que define o uso de terra comum.
#define COMEDI_REF 0
//! Parâmetro que define a frequência do sistema.
#define SYSTEM_FREQUENCY 60.00
//! Parâmetro que define quais harmônicos serão utilizados.
#define PRINT_N_HARMONICS 3
//! Parâmetro que define o tempo de atualização da visualização na tela.
#define REFRESH_TIME 1000 /*ms*/
//!Número máximo de oscilografias a serem armazenadas na Pillha
#define N_MAX_OSCILLOGRAM 20

#endif /* !PARAMETERS_H */

